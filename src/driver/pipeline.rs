//! The search of a composite with worker threads (see [`crate::Factorizer::threads`]).
//!
//! The sequential search is a sequence of steps: at each level, P-1 (or P+1) extended to the
//! bound of the level, then the curves of the level, one after the other, until a step finds a
//! factor. Here, the same steps (with the same curves: their parameters are drawn in the same
//! order) run on the workers of a [`Pool`] as soon as one is idle, ahead of the steps
//! reported: the curves of a level alongside its P-1 run and the last curves of the previous
//! level (P-1 runs resume one another: one at a time, in order). The factorization thread
//! reports the steps in their order (the events of the sequential search, but the durations),
//! and stops at the first step (in this order) that finds a factor, once all the steps before
//! it are reported: the steps after it, already running, are stopped and forgotten. The
//! result is the one of the sequential search.

use super::{
    Base2Form, Base2Mode, CurveOutcome, Engine, Error, Event, EventHandler, LEVELS, Method,
    PM1_B1_RATIO, PlusMinus, PlusMinusRun, Progress, Stage2Plan, TOP_LEVEL_ROUNDS, ecm_prob,
    pm1_b2, random_sigma, run_curve_timed, top_level,
};
use crate::{ecm::Param, parallel::Pool, stop::Stop};
use rug::{Integer, rand::RandState};
use std::{collections::VecDeque, sync::Arc, time::Duration};

/// How often the factorization thread checks the interruption flag of the
/// [`crate::Factorizer`] while it waits for the workers, to stop them.
const POLL: Duration = Duration::from_millis(1);

/// The parameter of a curve that panics (tests only).
#[cfg(test)]
pub(crate) const PANIC_SIGMA: u64 = 0x5eed_dead_beef;

/// The output of a job of the pool.
pub(crate) enum Output {
    /// A curve: its outcome, and the durations of its stages.
    Curve(CurveOutcome, [Duration; 2]),
    /// A run of P-1 (or P+1), with the state after it.
    PlusMinus(Box<(PlusMinus, PlusMinusRun)>),
}

/// What the curves of a level share.
struct Curves {
    n: Arc<Integer>,
    param: Param,
    /// Stage 1 multiplier, stage 2 plan and forms of the special reduction in both stages.
    k: Arc<Integer>,
    plan: Arc<Stage2Plan>,
    base2: [Option<Base2Form>; 2],
}

/// A step of the search.
enum Step {
    /// P-1 (or P+1) extended to `b1`.
    PlusMinus { b1: usize },
    /// The start of the curves of a level, reported by an [`Event::Level`].
    Level {
        digits: Option<u32>,
        b1: usize,
        b2: usize,
        curves: Option<usize>,
        done: usize,
    },
    /// A curve.
    Curve {
        curves: Arc<Curves>,
        index: usize,
        sigma: Integer,
    },
}

/// Where a step is.
enum State {
    /// To run (never for [`Step::Level`]).
    Waiting,
    /// Running on a worker.
    Running(usize),
    /// Ran (or nothing to run, for [`Step::Level`]).
    Done(Option<Output>),
    /// Not to run, nor reported: P-1 after it found all the factors at once.
    Skipped,
}

/// A step, with the progress of the search at this step.
struct Entry {
    step: Step,
    state: State,
    /// [`Progress::level`], [`Progress::rounds`] and [`Progress::pm_level`] at this step, and
    /// [`Progress::curves`] after it.
    level: usize,
    rounds: usize,
    pm_level: Option<usize>,
    curves: usize,
}

/// The steps of the sequential search, in order, from a [`Progress`] (with levels) or with fixed
/// bounds.
struct Steps<'r> {
    /// Fixed bounds and number of curves, if any.
    fixed: Option<(usize, usize, Option<usize>)>,
    top: usize,
    /// The progress of the sequential search (without P-1, see below).
    level: usize,
    rounds: usize,
    curves: usize,
    pm_level: Option<usize>,
    /// The bound of the last P-1 run, `None` without P-1 (it may still end when a run finds all
    /// the factors at once: its next steps are then skipped).
    pm_b1: Option<usize>,
    phase: Phase,
    /// Draws the parameters of the curves, in order.
    sigma: Sigma<'r>,
}

/// Where the sequence of steps is.
enum Phase {
    /// At the P-1 step of a level (maybe none).
    PlusMinus,
    /// At the start of the curves of a level.
    Level,
    /// At the curve `next` of a level (the last one is `last`, if any).
    Curves {
        curves: Arc<Curves>,
        next: usize,
        last: Option<usize>,
    },
    /// After the last step.
    End,
}

/// The parameters of the curves, drawn ahead of the curves reported: from copies of the random
/// state (or of the next parameter).
struct Sigma<'r> {
    rand: Option<RandState<'r>>,
    next: Option<Integer>,
}

impl Sigma<'_> {
    fn draw(&mut self, n: &Integer, param: Param) -> Integer {
        match (&mut self.next, &mut self.rand) {
            (Some(sigma), _) => {
                let next = Integer::from(&*sigma + 1u32);
                std::mem::replace(sigma, next)
            }
            (None, rand) => random_sigma(n, param, rand.as_mut().expect("a random state")),
        }
    }
}

impl<'r, H: EventHandler> Engine<'_, 'r, H> {
    /// Whether [`Engine::find`] runs the search of `n` with worker threads.
    pub(super) fn parallel(&self) -> bool {
        self.threads > 1
    }

    /// [`Engine::find`] with worker threads: with levels from `progress`, or with the fixed
    /// bounds and number of curves `fixed` (curves only).
    pub(super) fn find_parallel(
        &mut self,
        n: &Integer,
        progress: &mut Progress,
        base2: Option<Base2Form>,
        fixed: Option<(usize, usize, Option<usize>)>,
    ) -> Result<(Integer, Method), Error> {
        let top = top_level(n);
        let mut steps = Steps {
            fixed,
            top,
            level: progress.level,
            rounds: progress.rounds,
            curves: if fixed.is_some() { 0 } else { progress.curves },
            pm_level: progress.pm_level,
            pm_b1: progress
                .pm
                .as_ref()
                .map(PlusMinus::b1)
                .filter(|_| fixed.is_none()),
            phase: if fixed.is_some() {
                Phase::Level
            } else {
                Phase::PlusMinus
            },
            sigma: Sigma {
                // Only the generator of `rand_state` (a Factorizer's) can be copied.
                rand: self.sigma.is_none().then(|| self.rand.clone()),
                next: self.sigma.clone(),
            },
        };
        let n_shared = Arc::new(n.clone());
        let stop = self.events.stop;
        let poll = stop.polls_flag().then_some(POLL);
        let mut pool = self
            .pool
            .take()
            .unwrap_or_else(|| Pool::new(self.threads, stop.deadline()));
        let mut order: VecDeque<Entry> = VecDeque::new();
        // The P-1 state for the next P-1 step (`None` while one runs, `Some(None)` once P-1
        // found all the factors at once).
        let mut chain = Some(progress.pm.clone());
        // Position (from the front of `order`) of the first step that found a factor.
        let mut found: Option<usize> = None;
        let mut stopped = false;
        let mut reported_curves = 0;
        let result = 'search: loop {
            if !stopped && stop.requested() {
                stopped = true;
                pool.cancel_all();
            }
            // Start steps on the idle workers, in order.
            while let Some(worker) = pool.idle().filter(|_| !stopped) {
                let limit = found.unwrap_or(usize::MAX);
                let mut next = None;
                for (i, entry) in order.iter_mut().enumerate().take(limit) {
                    if !matches!(entry.state, State::Waiting) {
                        continue;
                    }
                    if let Step::PlusMinus { .. } = entry.step {
                        match &chain {
                            // The previous P-1 step still runs.
                            None => continue,
                            Some(None) => {
                                entry.state = State::Skipped;
                                continue;
                            }
                            Some(Some(_)) => {}
                        }
                    }
                    next = Some(i);
                    break;
                }
                let i = match next {
                    Some(i) => i,
                    None if found.is_some() => break,
                    None => match self.next_step(n, &n_shared, &mut steps, base2) {
                        Some(entry) => {
                            order.push_back(entry);
                            continue;
                        }
                        None => break,
                    },
                };
                let entry = &mut order[i];
                entry.state = State::Running(worker);
                match &entry.step {
                    Step::PlusMinus { b1 } => {
                        let mut pm = chain.take().flatten().expect("the P-1 state");
                        let (n, b1, max_memory) = (Arc::clone(&n_shared), *b1, self.max_memory);
                        let either = self.base2 == Base2Mode::Auto;
                        pool.submit(worker, move |stop: Stop<'_>| {
                            let plan = || {
                                let b2 = pm1_b2(&n, b1, max_memory, (base2, either));
                                Stage2Plan::cheapest(&n, (b1, b2), max_memory, base2, either)
                            };
                            let run =
                                PlusMinusRun::new(&mut pm, &n, (b1, plan), base2, H::ENABLED, stop);
                            Output::PlusMinus(Box::new((pm, run)))
                        });
                    }
                    Step::Curve { curves, sigma, .. } => {
                        let (curves, sigma) = (Arc::clone(curves), sigma.clone());
                        pool.submit(worker, move |stop: Stop<'_>| {
                            #[cfg(test)]
                            assert!(sigma != PANIC_SIGMA, "curve panicked");
                            let run = if H::ENABLED {
                                run_curve_timed::<true>
                            } else {
                                run_curve_timed::<false>
                            };
                            let c = &*curves;
                            let (outcome, times) =
                                run(&c.n, c.param, &sigma, &c.k, &c.plan, c.base2, stop);
                            Output::Curve(outcome, times)
                        });
                    }
                    Step::Level { .. } => unreachable!("nothing to run"),
                }
            }

            // Report the steps done, in order.
            while let Some(entry) = order.front() {
                match entry.state {
                    State::Waiting | State::Running(_) => break,
                    State::Skipped => {
                        order.pop_front();
                        found = found.map(|i| i - 1);
                        continue;
                    }
                    State::Done(_) => {}
                }
                let mut entry = order.pop_front().expect("a step");
                found = found.map(|i| i.saturating_sub(1));
                let State::Done(output) = std::mem::replace(&mut entry.state, State::Skipped)
                else {
                    unreachable!()
                };
                if stopped {
                    // After an interruption, only a factor found by the next step (which may
                    // have ended before it) is still reported.
                    let factor = match &output {
                        Some(Output::Curve(outcome, _)) => *outcome != CurveOutcome::Failed,
                        Some(Output::PlusMinus(run)) => run.1.found(n),
                        None => false,
                    };
                    if !factor {
                        self.events.interrupted = true;
                        break 'search Err(Error::Interrupted);
                    }
                }
                match (entry.step, output) {
                    (
                        Step::Level {
                            digits,
                            b1,
                            b2,
                            curves,
                            done,
                        },
                        _,
                    ) => {
                        let event = self.events.emit(Event::Level {
                            n,
                            digits,
                            b1,
                            b2,
                            curves,
                            done,
                        });
                        if let Err(error) = event {
                            break 'search Err(error);
                        }
                    }
                    (Step::Curve { index, sigma, .. }, Some(Output::Curve(outcome, times))) => {
                        reported_curves += 1;
                        let [stage1, stage2] = times;
                        let event = self.events.emit(Event::Curve {
                            n,
                            param: self.param,
                            sigma: &sigma,
                            index,
                            stage1,
                            stage2,
                        });
                        let found = match outcome {
                            CurveOutcome::Setup(g) => (g, Method::EcmSetup),
                            CurveOutcome::Stage1(g) => (g, Method::EcmStage1),
                            CurveOutcome::Stage2(g) => (g, Method::EcmStage2),
                            CurveOutcome::Failed => match event {
                                Ok(()) => continue,
                                Err(error) => break 'search Err(error),
                            },
                        };
                        (progress.level, progress.rounds) = (entry.level, entry.rounds);
                        (progress.pm_level, progress.curves) = (entry.pm_level, entry.curves);
                        break 'search Ok(found);
                    }
                    (Step::PlusMinus { b1 }, Some(Output::PlusMinus(run))) => {
                        let (pm, run) = *run;
                        let all = run.all(n);
                        let report = run.report(&pm, n, b1, &mut self.events);
                        progress.pm = (!all).then_some(pm);
                        match report {
                            Ok((None, _)) => {}
                            Ok((Some(found), _)) => {
                                (progress.level, progress.rounds) = (entry.level, entry.rounds);
                                (progress.pm_level, progress.curves) =
                                    (entry.pm_level, entry.curves);
                                break 'search Ok(found);
                            }
                            Err(error) => break 'search Err(error),
                        }
                    }
                    _ => unreachable!("the output of another step"),
                }
            }

            if !pool.busy() {
                // Nothing runs: all the steps ran (or the search is interrupted).
                debug_assert!(order.is_empty() || stopped);
                break Err(if self.events.interrupted() {
                    Error::Interrupted
                } else {
                    Error::ECMFailed
                });
            }
            let Some((worker, output)) = pool.recv(poll) else {
                continue;
            };
            let Some(i) = order
                .iter()
                .position(|entry| matches!(entry.state, State::Running(w) if w == worker))
            else {
                unreachable!("a step runs on the worker")
            };
            let factor = match &output {
                Output::Curve(outcome, _) => *outcome != CurveOutcome::Failed,
                Output::PlusMinus(run) => {
                    let (pm, run) = &**run;
                    chain = Some((!run.all(n)).then(|| pm.clone()));
                    run.found(n)
                }
            };
            order[i].state = State::Done(Some(output));
            if factor && found.is_none_or(|found| i < found) {
                found = Some(i);
                // The steps after it are useless.
                for entry in order.iter().skip(i + 1) {
                    if let State::Running(worker) = entry.state {
                        pool.cancel(worker);
                    }
                }
            }
        };
        pool.stop_all();
        self.pool = Some(pool);

        // The random state (or the next parameter) as after the curves reported.
        match &mut self.sigma {
            Some(sigma) => *sigma += reported_curves,
            None => {
                for _ in 0..reported_curves {
                    random_sigma(n, self.param, self.rand);
                }
            }
        }
        result
    }

    /// The next step of the sequential search, if any.
    fn next_step(
        &mut self,
        n: &Integer,
        n_shared: &Arc<Integer>,
        steps: &mut Steps<'r>,
        base2: Option<Base2Form>,
    ) -> Option<Entry> {
        loop {
            let index = steps.level.min(steps.top);
            match &mut steps.phase {
                Phase::PlusMinus => {
                    steps.phase = Phase::Level;
                    if steps.pm_level >= Some(index) {
                        continue;
                    }
                    steps.pm_level = Some(index);
                    let b1 = LEVELS[index].b1 * PM1_B1_RATIO;
                    if steps.pm_b1.is_some_and(|pm_b1| b1 > pm_b1) {
                        steps.pm_b1 = Some(b1);
                        let step = Step::PlusMinus { b1 };
                        return Some(steps.entry(step, State::Waiting, steps.curves));
                    }
                }
                Phase::Level => {
                    let (b1, b2, digits, curves) = match steps.fixed {
                        Some((b1, b2, curves)) => (b1, b2, None, curves),
                        None => {
                            let level = &LEVELS[index];
                            let (plan, _) = self.plan(n, level.b1, level.b2, base2);
                            let prob =
                                ecm_prob(level.b1 as f64, plan.b2() as f64, level.digits.into());
                            let curves = (1.0 / prob).ceil().max(1.0) as usize;
                            (level.b1, level.b2, Some(level.digits), Some(curves))
                        }
                    };
                    let done = steps.curves;
                    if curves.is_some_and(|curves| done >= curves) {
                        // Resumed on a cofactor after all the curves of the level.
                        steps.end_level();
                        continue;
                    }
                    let k = self.multiplier(b1);
                    let (plan, base2_stage2) = self.plan(n, b1, b2, base2);
                    let b2 = plan.b2();
                    steps.phase = Phase::Curves {
                        curves: Arc::new(Curves {
                            n: Arc::clone(n_shared),
                            param: self.param,
                            k,
                            plan,
                            base2: [base2, base2_stage2],
                        }),
                        next: done + 1,
                        last: curves,
                    };
                    let step = Step::Level {
                        digits,
                        b1,
                        b2,
                        curves,
                        done,
                    };
                    return Some(steps.entry(step, State::Done(None), done));
                }
                Phase::Curves { curves, next, last } => {
                    if last.is_some_and(|last| *next > last) {
                        steps.end_level();
                        continue;
                    }
                    let (curves, index) = (Arc::clone(curves), *next);
                    *next += 1;
                    let sigma = steps.sigma.draw(n, self.param);
                    let step = Step::Curve {
                        curves,
                        index,
                        sigma,
                    };
                    return Some(steps.entry(step, State::Waiting, index));
                }
                Phase::End => return None,
            }
        }
    }
}

impl Steps<'_> {
    /// A step at the current progress, with `curves` curves of the level run after it.
    fn entry(&self, step: Step, state: State, curves: usize) -> Entry {
        Entry {
            step,
            state,
            level: self.level,
            rounds: self.rounds,
            pm_level: self.pm_level,
            curves,
        }
    }

    /// After the curves of a level: the next level (or the next round of the last one).
    fn end_level(&mut self) {
        self.curves = 0;
        let index = self.level.min(self.top);
        self.phase = if self.fixed.is_some() {
            Phase::End
        } else if index < self.top {
            self.level = index + 1;
            Phase::PlusMinus
        } else {
            self.rounds += 1;
            if self.rounds >= TOP_LEVEL_ROUNDS {
                Phase::End
            } else {
                Phase::PlusMinus
            }
        };
    }
}
