//! Curves run in parallel (see [`crate::Factorizer::threads`]): a pool of worker threads, each
//! running the curves the factorization thread gives it, one at a time.
//!
//! The factorization thread draws the parameters of the curves in their order, hands them to
//! the idle workers, and receives the outcomes: it reports them in the order of the curves, and
//! stops at the first curve (in that order) that finds a factor, once all the curves before it
//! ran. The results are then those of the same curves run one after the other.

use crate::{
    base2::Base2Form,
    ecm::{CurveOutcome, Param, run_curve_timed},
    stage2::Stage2Plan,
    stop::Stop,
};
use rug::Integer;
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

/// What the curves of a level share: the number, the stage 1 multiplier, the stage 2 plan and
/// the forms of the special reduction in both stages.
pub(crate) struct Curves {
    pub(crate) n: Arc<Integer>,
    pub(crate) param: Param,
    pub(crate) k: Arc<Integer>,
    pub(crate) plan: Arc<Stage2Plan>,
    pub(crate) base2: [Option<Base2Form>; 2],
    /// Whether the durations of the stages are measured.
    pub(crate) timed: bool,
    /// When the curves must stop.
    pub(crate) deadline: Option<Instant>,
}

/// A curve to run.
struct Job {
    curves: Arc<Curves>,
    index: usize,
    sigma: Integer,
}

/// The outcome of a curve, with the durations of its stages.
pub(crate) type Timed = (CurveOutcome, [Duration; 2]);

/// The outcome of a curve, or the panic of its worker.
type Outcome = thread::Result<Timed>;

/// A curve run by a worker.
struct Done {
    worker: usize,
    index: usize,
    outcome: Outcome,
}

/// A worker thread, with the curve it runs.
struct Worker {
    jobs: Option<Sender<Job>>,
    /// Stops its curve (checked as the interruption flag of [`crate::Factorizer`] is).
    cancel: Arc<AtomicBool>,
    /// The index and parameter of its curve, if it runs one.
    busy: Option<(usize, Integer)>,
    handle: Option<JoinHandle<()>>,
}

/// Worker threads (started at the first use), which run curves.
pub(crate) struct Pool {
    threads: usize,
    workers: Vec<Worker>,
    done_tx: Sender<Done>,
    done_rx: Receiver<Done>,
}

impl Pool {
    /// A pool of `threads` workers, not started yet.
    pub(crate) fn new(threads: usize) -> Self {
        let (done_tx, done_rx) = mpsc::channel();
        Self {
            threads,
            workers: Vec::new(),
            done_tx,
            done_rx,
        }
    }

    fn start(&mut self) {
        while self.workers.len() < self.threads {
            let id = self.workers.len();
            let (jobs, rx) = mpsc::channel::<Job>();
            let cancel = Arc::new(AtomicBool::new(false));
            let done = self.done_tx.clone();
            let flag = Arc::clone(&cancel);
            let handle = thread::Builder::new()
                .name(format!("ecm-curves-{id}"))
                .spawn(move || {
                    for job in rx {
                        let curves = &*job.curves;
                        let stop = Stop::new(Some(&flag), curves.deadline);
                        let outcome = panic::catch_unwind(AssertUnwindSafe(|| {
                            #[cfg(test)]
                            assert!(job.sigma != tests::PANIC_SIGMA, "curve panicked");
                            let run = if curves.timed {
                                run_curve_timed::<true>
                            } else {
                                run_curve_timed::<false>
                            };
                            let (n, k, plan) = (&curves.n, &curves.k, &curves.plan);
                            run(n, curves.param, &job.sigma, k, plan, curves.base2, stop)
                        }));
                        let done_msg = Done {
                            worker: id,
                            index: job.index,
                            outcome,
                        };
                        if done.send(done_msg).is_err() {
                            break;
                        }
                    }
                })
                .expect("failed to start a thread");
            self.workers.push(Worker {
                jobs: Some(jobs),
                cancel,
                busy: None,
                handle: Some(handle),
            });
        }
    }

    /// Runs the curves `from, from + 1, ...` (up to `last` if any) of `curves`, with the
    /// parameters given by `sigma` (called in the order of the curves), until `report` (called
    /// in the order of the curves with each outcome) returns `true`, or `stop` is requested.
    ///
    /// Returns the number of curves reported. A curve is reported only once all the curves
    /// before it were; after a stop, only the next curve may still be, if it found a factor.
    /// All the curves running are stopped before returning.
    pub(crate) fn run(
        &mut self,
        curves: &Arc<Curves>,
        (from, last): (usize, Option<usize>),
        mut sigma: impl FnMut() -> Integer,
        stop: Stop<'_>,
        mut report: impl FnMut(usize, Integer, Timed) -> bool,
    ) -> usize {
        self.start();
        // Outcomes received but not reported yet (the curves before them still run), by index.
        let mut received: Vec<Option<(Integer, Timed)>> = Vec::new();
        let mut next_index = from;
        let mut next_report = from;
        // Lowest index of a curve that found a factor (the curves after it are useless).
        let mut found = usize::MAX;
        let mut stopped = false;
        let mut finished = false;
        let poll = stop.polls_flag().then_some(POLL);
        loop {
            if !stopped && !finished {
                for worker in &mut self.workers {
                    if worker.busy.is_some() {
                        continue;
                    }
                    if next_index >= found || last.is_some_and(|last| next_index > last) {
                        break;
                    }
                    let job_sigma = sigma();
                    worker.cancel.store(false, Ordering::Relaxed);
                    let job = Job {
                        curves: Arc::clone(curves),
                        index: next_index,
                        sigma: job_sigma.clone(),
                    };
                    worker.busy = Some((next_index, job_sigma));
                    // The worker only stops when the pool is dropped.
                    let _ = worker.jobs.as_ref().map(|jobs| jobs.send(job));
                    next_index += 1;
                }
            }
            if self.workers.iter().all(|worker| worker.busy.is_none()) {
                break;
            }
            let done = match poll {
                Some(poll) => match self.done_rx.recv_timeout(poll) {
                    Ok(done) => Some(done),
                    Err(RecvTimeoutError::Timeout) => None,
                    Err(RecvTimeoutError::Disconnected) => unreachable!("the pool keeps a sender"),
                },
                None => Some(self.done_rx.recv().expect("the pool keeps a sender")),
            };
            if !stopped && stop.requested() {
                stopped = true;
                self.cancel_all();
            }
            let Some(Done {
                worker,
                index,
                outcome,
            }) = done
            else {
                continue;
            };
            let (_, sigma) = self.workers[worker].busy.take().expect("a worker is busy");
            let outcome = match outcome {
                Ok(outcome) => outcome,
                Err(payload) => {
                    self.cancel_all();
                    self.wait_idle();
                    panic::resume_unwind(payload);
                }
            };
            if finished || index >= found {
                continue;
            }
            if outcome.0 == CurveOutcome::Failed {
                if stopped && index == next_report {
                    // It may have stopped early: not reported, and nothing is after it.
                    finished = true;
                    self.cancel_all();
                    continue;
                }
            } else {
                found = index;
                for worker in &self.workers {
                    if worker.busy.as_ref().is_some_and(|(i, _)| *i > index) {
                        worker.cancel.store(true, Ordering::Relaxed);
                    }
                }
            }
            let slot = index - next_report;
            if received.len() <= slot {
                received.resize_with(slot + 1, || None);
            }
            received[slot] = Some((sigma, outcome));
            while let Some(Some(_)) = received.first() {
                let (sigma, outcome) = received.remove(0).expect("received");
                let index = next_report;
                next_report += 1;
                if report(index, sigma, outcome) || index == found {
                    finished = true;
                    self.cancel_all();
                    break;
                }
            }
        }
        next_report - from
    }

    fn cancel_all(&self) {
        for worker in &self.workers {
            worker.cancel.store(true, Ordering::Relaxed);
        }
    }

    /// Waits for the curves running to stop (they must have been cancelled).
    fn wait_idle(&mut self) {
        while self.workers.iter().any(|worker| worker.busy.is_some()) {
            let done = self.done_rx.recv().expect("the pool keeps a sender");
            self.workers[done.worker].busy = None;
        }
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        self.cancel_all();
        for worker in &mut self.workers {
            worker.jobs = None;
        }
        for worker in &mut self.workers {
            if let Some(handle) = worker.handle.take() {
                // The workers catch the panics of the curves.
                let _ = handle.join();
            }
        }
    }
}

/// How often the factorization thread checks the interruption flag of the
/// [`crate::Factorizer`] while the workers run curves, to stop them.
const POLL: Duration = Duration::from_millis(1);

#[cfg(test)]
mod tests {
    use crate::{Error, Factorizer};
    use rug::Integer;
    use std::{
        panic,
        time::{Duration, Instant},
    };

    /// The parameter of a curve that panics.
    pub(super) const PANIC_SIGMA: u64 = 0x5eed_dead_beef;

    #[test]
    fn panics_are_propagated() {
        // 20-digit factors: many curves, so the one that panics is not the first.
        let p = Integer::from(10u64.pow(19)).next_prime();
        let n = Integer::from(&p * &p).next_prime() * p;
        for threads in [2, 4] {
            let start = Instant::now();
            let result = panic::catch_unwind(|| {
                Factorizer::new()
                    .b1(11_000)
                    .sigma(Integer::from(PANIC_SIGMA - 5))
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
            .sigma(Integer::from(PANIC_SIGMA - 5))
            .threads(4)
            .factor(&n);
        assert_eq!(result, Err(Error::ECMFailed));
    }
}
