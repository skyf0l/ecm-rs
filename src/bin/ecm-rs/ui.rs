//! Diagnostics on stderr, built from the events of the factorization: `-v` lines (worded as
//! GMP-ECM's) and the progress bar.

use ecm::{Algorithm, Event, Method, Param};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rug::Integer;
use std::{fmt::Write as _, time::Duration};

/// Number of decimal digits of `n > 0`.
pub fn digits(n: &Integer) -> usize {
    n.to_string_radix(10).trim_start_matches('-').len()
}

/// `11000` as `11k`, `1873422` as `1.9M`.
pub fn short(x: usize) -> String {
    let x = x as f64;
    for (unit, scale) in [("T", 1e12), ("G", 1e9), ("M", 1e6), ("k", 1e3)] {
        if x >= scale {
            let v = x / scale;
            return if v >= 100.0 || (v - v.round()).abs() < 0.05 {
                format!("{v:.0}{unit}")
            } else {
                format!("{v:.1}{unit}")
            };
        }
    }
    format!("{x}")
}

/// A duration as GMP-ECM prints it (`12ms`), with a decimal below 10 ms (`0.4ms`).
fn ms(d: Duration) -> String {
    let ms = d.as_secs_f64() * 1e3;
    if ms < 10.0 {
        format!("{ms:.1}ms")
    } else {
        format!("{ms:.0}ms")
    }
}

/// A duration for the progress bar: `850ms`, `12s`, `3m20s`, `2h05m`.
fn human(d: Duration) -> String {
    let s = d.as_secs_f64();
    if s < 1.0 {
        format!("{}ms", d.as_millis())
    } else if s < 60.0 {
        format!("{s:.0}s")
    } else if s < 3600.0 {
        format!("{}m{:02}s", (s / 60.0) as u64, s as u64 % 60)
    } else {
        format!("{}h{:02}m", (s / 3600.0) as u64, (s as u64 / 60) % 60)
    }
}

/// The level being run.
#[derive(Default)]
struct Level {
    /// Whether the bounds are those of a factor size (not fixed).
    by_size: bool,
    b1: usize,
    b2: usize,
    curves: Option<usize>,
    /// Curves timed at this level, and their total duration.
    timed: u32,
    time: Duration,
}

/// Reports the events of one factorization on stderr.
pub struct Ui {
    verbose: u8,
    /// P-1 or P+1, if it runs first (not with fixed bounds, but when it runs alone).
    first: Option<Algorithm>,
    /// The seed of P-1 or P+1, if not the default one.
    x0: Option<(Integer, Integer)>,
    bar: Option<ProgressBar>,
    level: Level,
    /// The last curve run: `(param, sigma)`, for the factor it may have found.
    last_curve: Option<(Param, Integer)>,
}

impl Ui {
    /// `verbose` lines (0: none), and a progress bar if `progress`; `first`: P-1 or P+1 if it
    /// runs first; `x0`: its seed, if given.
    pub fn new(
        verbose: u8,
        progress: bool,
        first: Option<Algorithm>,
        x0: Option<(Integer, Integer)>,
    ) -> Self {
        let bar = progress.then(|| {
            let bar = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
            bar.set_style(spinner_style());
            bar.set_message("trial division");
            bar.enable_steady_tick(Duration::from_millis(120));
            bar
        });
        Self {
            verbose,
            first,
            x0,
            bar,
            level: Level::default(),
            last_curve: None,
        }
    }

    /// Prints `line` on stderr, above the progress bar.
    pub fn line(&self, line: &str) {
        match &self.bar {
            Some(bar) => bar.println(line),
            None => eprintln!("{line}"),
        }
    }

    /// Whether the factors found are reported.
    fn show_factors(&self) -> bool {
        self.verbose > 0 || self.bar.is_some()
    }

    /// Removes the progress bar.
    pub fn finish(&self) {
        if let Some(bar) = &self.bar {
            bar.finish_and_clear();
        }
    }

    pub fn handle(&mut self, event: &Event<'_>) {
        match *event {
            Event::TrialDivision {
                factors, cofactor, ..
            } => {
                if self.verbose > 0 {
                    let mut primes: Vec<_> = factors.iter().collect();
                    primes.sort();
                    let found = if primes.is_empty() {
                        "no factor".to_string()
                    } else {
                        product(primes.iter().map(|&(p, &e)| (p, e)))
                    };
                    let rest = if *cofactor == 1 {
                        String::new()
                    } else {
                        format!(", cofactor has {} digits", digits(cofactor))
                    };
                    self.line(&format!("Trial division below 2^16: {found}{rest}"));
                }
                let next = self.first.unwrap_or(Algorithm::Ecm);
                self.spinner(format!("{next} on C{}", digits(cofactor)));
            }
            Event::Base2 { k, .. } => {
                if self.verbose > 0 {
                    let sign = if k > 0 { '+' } else { '-' };
                    let k = k.unsigned_abs();
                    self.line(&format!(
                        "Using special division for factor of 2^{k}{sign}1"
                    ));
                }
            }
            Event::Pm1 {
                n,
                b1,
                b2,
                stage1,
                stage2,
                ..
            }
            | Event::Pp1 {
                n,
                b1,
                b2,
                stage1,
                stage2,
                ..
            } => {
                if self.verbose > 0 {
                    let (name, default) = match event {
                        Event::Pp1 { .. } => ("P+1", "2/7"),
                        _ => ("P-1", "3"),
                    };
                    let mut line = format!("{name} on C{}: B1={b1}", digits(n));
                    if let Some(b2) = b2 {
                        let _ = write!(line, ", B2={b2}");
                    }
                    let _ = match &self.x0 {
                        Some((num, den)) if *den == 1 => write!(line, ", x0={num}"),
                        Some((num, den)) => write!(line, ", x0={num}/{den}"),
                        None => write!(line, ", x0={default}"),
                    };
                    let _ = write!(line, ": Step 1 took {}", ms(stage1));
                    if b2.is_some() {
                        let _ = write!(line, ", Step 2 took {}", ms(stage2));
                    }
                    self.line(&line);
                }
            }
            Event::Level {
                n,
                digits: level_digits,
                b1,
                b2,
                curves,
                done,
                ..
            } => {
                self.level = Level {
                    by_size: level_digits.is_some(),
                    b1,
                    b2,
                    curves,
                    ..Level::default()
                };
                if self.verbose > 0 {
                    let mut line = format!("Using B1={b1}, B2={b2} on C{}", digits(n));
                    if let Some(d) = level_digits {
                        let _ = write!(line, ", factors of {d} digits");
                    }
                    if let Some(curves) = curves {
                        let what = if level_digits.is_some() {
                            "expected curves"
                        } else {
                            "curves"
                        };
                        let _ = write!(line, ": {what} {curves}");
                        if done > 0 {
                            let _ = write!(line, " ({done} already run)");
                        }
                    }
                    self.line(&line);
                }
                if let Some(bar) = &self.bar {
                    let name = match level_digits {
                        Some(d) => format!("ECM p{d}"),
                        None => "ECM".to_string(),
                    };
                    bar.set_prefix(format!("{name} B1={} B2={}", short(b1), short(b2)));
                    bar.set_message("");
                    match curves {
                        Some(curves) => {
                            bar.set_style(bar_style());
                            bar.set_length(curves as u64);
                        }
                        None => bar.set_style(count_style()),
                    }
                    bar.reset();
                    bar.set_position(done as u64);
                }
            }
            Event::Curve {
                param,
                sigma,
                index,
                stage1,
                stage2,
                ..
            } => {
                self.last_curve = Some((param, sigma.clone()));
                if self.verbose > 0 {
                    let of = self
                        .level
                        .curves
                        .map_or(String::new(), |curves| format!("/{curves}"));
                    let mut line = format!(
                        "Curve {index}{of}: sigma={param}:{sigma}, Step 1 took {}",
                        ms(stage1)
                    );
                    if stage2 > Duration::ZERO {
                        let _ = write!(line, ", Step 2 took {}", ms(stage2));
                    }
                    self.line(&line);
                }
                let level = &mut self.level;
                level.timed += 1;
                level.time += stage1 + stage2;
                if let Some(bar) = &self.bar {
                    bar.set_position(index as u64);
                    match level.curves {
                        Some(curves) if index >= curves && level.by_size => {
                            // P-1 comes before the next level.
                            self.spinner("P-1".to_string());
                        }
                        Some(curves) => {
                            let mean = level.time / level.timed;
                            let expected = mean.saturating_mul(curves as u32);
                            let expected = human(expected);
                            bar.set_message(format!("(expected ~{expected} at this level)"));
                        }
                        None => {}
                    }
                }
            }
            Event::Factor {
                n, factor, method, ..
            } => {
                if !self.show_factors() {
                    return;
                }
                let how = match (method, &self.last_curve) {
                    (Method::EcmSetup | Method::EcmStage1 | Method::EcmStage2, Some((p, s))) => {
                        let bounds = format!("B1={}, B2={}", self.level.b1, self.level.b2);
                        format!("{method} ({bounds}, sigma={p}:{s})")
                    }
                    _ => method.to_string(),
                };
                self.line(&format!("********** Factor found by {how}: {factor}"));
                if matches!(method, Method::Pp1Stage1 | Method::Pp1Stage2)
                    && self.verbose > 0
                    && is_prime(factor)
                {
                    // As GMP-ECM: x0^2 - 4 is a square modulo p, the group has order p - 1
                    // (unless p divides the denominator: x0 is not defined modulo p).
                    let (num, den) = self.x0.clone().unwrap_or((2.into(), 7.into()));
                    let d = Integer::from(&num * &num) - Integer::from(&den * &den) * 4u32;
                    if !den.is_divisible(factor) && d.jacobi(factor) == 1 {
                        self.line("[factor found by P-1]");
                    }
                }
                if self.verbose > 0 {
                    let cofactor = Integer::from(n / factor);
                    let kind = |x: &Integer| {
                        if is_prime(x) { "prime" } else { "composite" }
                    };
                    let (fk, ck) = (kind(factor), kind(&cofactor));
                    self.line(&format!(
                        "Found {fk} factor of {} digits: {factor}",
                        digits(factor)
                    ));
                    let cap = |s: &str| format!("{}{}", s[..1].to_uppercase(), &s[1..]);
                    self.line(&format!(
                        "{} cofactor {cofactor} has {} digits",
                        cap(ck),
                        digits(&cofactor)
                    ));
                }
            }
            Event::Prime { p, exponent, .. } if self.verbose > 1 => {
                self.line(&format!("Prime factor: {}", product([(p, exponent)])));
            }
            _ => {}
        }
    }

    /// Shows the spinner with `msg`, between levels.
    fn spinner(&self, msg: String) {
        if let Some(bar) = &self.bar {
            bar.set_message(msg);
            bar.set_style(spinner_style());
        }
    }
}

/// Whether `n` is (probably) prime, as the library tests it.
pub fn is_prime(n: &Integer) -> bool {
    n.is_probably_prime(25) != rug::integer::IsPrime::No
}

/// `p1^e1 * p2 * ...`.
pub fn product<'a>(factors: impl IntoIterator<Item = (&'a Integer, usize)>) -> String {
    let mut s = String::new();
    for (p, e) in factors {
        if !s.is_empty() {
            s.push_str(" * ");
        }
        let _ = write!(s, "{p}");
        if e > 1 {
            let _ = write!(s, "^{e}");
        }
    }
    s
}

fn spinner_style() -> ProgressStyle {
    ProgressStyle::with_template("{spinner} {msg} {elapsed}").expect("valid template")
}

fn bar_style() -> ProgressStyle {
    ProgressStyle::with_template("{prefix} [{wide_bar}] {pos}/{len} curves {elapsed} {msg}")
        .expect("valid template")
        .progress_chars("=> ")
}

fn count_style() -> ProgressStyle {
    ProgressStyle::with_template("{spinner} {prefix} {pos} curves {elapsed}")
        .expect("valid template")
}
