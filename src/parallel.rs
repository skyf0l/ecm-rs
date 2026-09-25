//! Worker threads for [`crate::Factorizer::threads`]: a pool running jobs (curves, P-1 or P+1)
//! given by the factorization thread, which also receives their outputs (see the parallel search
//! in [`crate::driver`]).

use crate::stop::Stop;
use std::{
    panic::{self, AssertUnwindSafe},
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
        mpsc::{self, Receiver, RecvTimeoutError, Sender},
    },
    thread::{self, JoinHandle},
    time::{Duration, Instant},
};

/// A job: a computation stopping early when its [`Stop`] is requested.
type Job<T> = Box<dyn FnOnce(Stop<'_>) -> T + Send>;

/// The output of a job (or the panic of its worker), by the worker that ran it.
type Done<T> = (usize, thread::Result<T>);

/// A worker thread.
struct Worker<T> {
    jobs: Option<Sender<Job<T>>>,
    /// Stops its job (checked as the interruption flag of [`crate::Factorizer`] is).
    cancel: Arc<AtomicBool>,
    busy: bool,
    handle: Option<JoinHandle<()>>,
}

/// Worker threads (started at the first use) running jobs with the output `T`, one at a time
/// each. Dropping the pool stops the jobs running, and waits for the threads to end.
pub(crate) struct Pool<T> {
    threads: usize,
    /// When the jobs must stop, as the factorization.
    deadline: Option<Instant>,
    workers: Vec<Worker<T>>,
    done_tx: Sender<Done<T>>,
    done_rx: Receiver<Done<T>>,
}

impl<T: Send + 'static> Pool<T> {
    /// A pool of `threads` workers (not started yet) whose jobs stop at `deadline`.
    pub(crate) fn new(threads: usize, deadline: Option<Instant>) -> Self {
        let (done_tx, done_rx) = mpsc::channel();
        Self {
            threads,
            deadline,
            workers: Vec::new(),
            done_tx,
            done_rx,
        }
    }

    fn start(&mut self) {
        while self.workers.len() < self.threads {
            let id = self.workers.len();
            let (jobs, rx) = mpsc::channel::<Job<T>>();
            let cancel = Arc::new(AtomicBool::new(false));
            let (done, flag, deadline) = (self.done_tx.clone(), Arc::clone(&cancel), self.deadline);
            let handle = thread::Builder::new()
                .name(format!("ecm-{id}"))
                .spawn(move || {
                    for job in rx {
                        let stop = Stop::new(Some(&flag), deadline);
                        let output = panic::catch_unwind(AssertUnwindSafe(|| job(stop)));
                        if done.send((id, output)).is_err() {
                            break;
                        }
                    }
                })
                .expect("failed to start a thread");
            self.workers.push(Worker {
                jobs: Some(jobs),
                cancel,
                busy: false,
                handle: Some(handle),
            });
        }
    }

    /// An idle worker, if any (starting the workers at the first call).
    pub(crate) fn idle(&mut self) -> Option<usize> {
        self.start();
        self.workers.iter().position(|worker| !worker.busy)
    }

    /// Whether a job runs.
    pub(crate) fn busy(&self) -> bool {
        self.workers.iter().any(|worker| worker.busy)
    }

    /// Runs `job` on the idle `worker`.
    pub(crate) fn submit(
        &mut self,
        worker: usize,
        job: impl FnOnce(Stop<'_>) -> T + Send + 'static,
    ) {
        let worker = &mut self.workers[worker];
        debug_assert!(!worker.busy);
        worker.cancel.store(false, Ordering::Relaxed);
        worker.busy = true;
        if let Some(jobs) = &worker.jobs {
            // The workers only end when the pool is dropped.
            let _ = jobs.send(Box::new(job));
        }
    }

    /// Stops the job of `worker` (if it still runs, its output is received anyway).
    pub(crate) fn cancel(&self, worker: usize) {
        self.workers[worker].cancel.store(true, Ordering::Relaxed);
    }

    /// Stops all the jobs.
    pub(crate) fn cancel_all(&self) {
        for worker in &self.workers {
            worker.cancel.store(true, Ordering::Relaxed);
        }
    }

    /// The output of the next job to end, with its worker (now idle), waiting at most `timeout`
    /// if any. Re-raises the panic of a job, once the other jobs are stopped.
    ///
    /// # Panics
    ///
    /// If no job runs, or a job panicked.
    pub(crate) fn recv(&mut self, timeout: Option<Duration>) -> Option<(usize, T)> {
        assert!(self.busy(), "no job runs");
        let (worker, output) = match timeout {
            Some(timeout) => match self.done_rx.recv_timeout(timeout) {
                Ok(done) => done,
                Err(RecvTimeoutError::Timeout) => return None,
                Err(RecvTimeoutError::Disconnected) => unreachable!("the pool keeps a sender"),
            },
            None => self.done_rx.recv().expect("the pool keeps a sender"),
        };
        self.workers[worker].busy = false;
        match output {
            Ok(output) => Some((worker, output)),
            Err(payload) => {
                self.stop_all();
                panic::resume_unwind(payload);
            }
        }
    }

    /// Stops the jobs running, and waits for them to end (discarding their outputs).
    pub(crate) fn stop_all(&mut self) {
        self.cancel_all();
        while self.busy() {
            let (worker, _) = self.done_rx.recv().expect("the pool keeps a sender");
            self.workers[worker].busy = false;
        }
    }
}

impl<T> Drop for Pool<T> {
    fn drop(&mut self) {
        for worker in &mut self.workers {
            worker.cancel.store(true, Ordering::Relaxed);
            worker.jobs = None;
        }
        for worker in &mut self.workers {
            if let Some(handle) = worker.handle.take() {
                // The workers catch the panics of the jobs.
                let _ = handle.join();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Error, Factorizer};
    use rug::Integer;

    #[test]
    fn pool() {
        let mut pool = Pool::<usize>::new(3, None);
        for i in 0..3 {
            let worker = pool.idle().unwrap();
            pool.submit(worker, move |stop: Stop<'_>| {
                while !stop.requested() {
                    thread::sleep(Duration::from_millis(1));
                }
                i
            });
        }
        assert_eq!(pool.idle(), None);
        pool.cancel(1);
        assert_eq!(pool.recv(None), Some((1, 1)));
        assert_eq!(pool.recv(Some(Duration::from_millis(5))), None);
        pool.stop_all();
        assert!(!pool.busy());
    }

    #[test]
    fn panics_are_propagated() {
        // 20-digit factors: many curves, so the one that panics is not the first.
        let p = Integer::from(10u64.pow(19)).next_prime();
        let n = Integer::from(&p * &p).next_prime() * p;
        let sigma = Integer::from(crate::driver::PANIC_SIGMA - 5);
        for threads in [2, 4] {
            let start = Instant::now();
            let result = panic::catch_unwind(|| {
                Factorizer::new()
                    .b1(11_000)
                    .sigma(sigma.clone())
                    .threads(threads)
                    .factor(&n)
            });
            let payload = result.expect_err("the panic of a curve is propagated");
            assert_eq!(payload.downcast_ref::<&str>(), Some(&"curve panicked"));
            assert!(start.elapsed() < Duration::from_secs(5));
        }
        // Before the curve that panics: no panic.
        let result = Factorizer::new()
            .b1(11_000)
            .curves(5)
            .sigma(sigma)
            .threads(4)
            .factor(&n);
        assert_eq!(result, Err(Error::ECMFailed));
    }
}
