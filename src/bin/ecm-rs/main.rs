//! `ecm-rs`: factors integers with the `ecm` crate, with GMP-ECM-like options.

mod args;
// The library's build choices, for `--printconfig`.
#[path = "../../config.rs"]
mod config;
mod expr;
mod ui;

use args::{Cli, status};
use clap::{Parser, error::ErrorKind};
use ecm::{Algorithm, Error, Factorizer, Param};
use rug::Integer;
use serde_json::{Value, json};
use std::{
    borrow::Cow,
    collections::HashMap,
    io::{BufRead, IsTerminal, Write},
    ops::ControlFlow,
    process::ExitCode,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    time::Instant,
};
use ui::{Ui, digits, is_prime, product};

fn main() -> ExitCode {
    let cli = match Cli::try_parse() {
        Ok(cli) => cli,
        Err(e) => {
            let _ = e.print();
            return ExitCode::from(match e.kind() {
                ErrorKind::DisplayHelp | ErrorKind::DisplayVersion => 0,
                _ => status::USAGE,
            });
        }
    };
    if cli.printconfig {
        print_config();
        return ExitCode::SUCCESS;
    }
    let factorizer = match factorizer(&cli) {
        Ok(factorizer) => factorizer,
        Err(msg) => {
            eprintln!("ecm-rs: {msg}");
            return ExitCode::from(status::USAGE);
        }
    };

    let interrupted = Arc::new(AtomicBool::new(false));
    let busy = Arc::new(AtomicBool::new(false));
    {
        let (flag, busy) = (interrupted.clone(), busy.clone());
        // Ctrl-C interrupts the factorization in progress (its partial result is printed); at
        // any other time (waiting for the standard input, testing a prime), or a second time,
        // it exits at once. No handler (an error): Ctrl-C kills the process.
        let _ = ctrlc::set_handler(move || {
            if !busy.load(Ordering::SeqCst) || flag.swap(true, Ordering::SeqCst) {
                std::process::exit(status::INTERRUPTED.into());
            }
        });
    }
    let factorizer = factorizer.interrupt_flag(interrupted.clone());

    let mut run = Run {
        cli: &cli,
        factorizer,
        interrupted,
        busy,
        deadline: cli.total_timeout.map(|t| Instant::now() + t),
        progress: !cli.quiet
            && !cli.json
            && !cli.no_progress
            && !cli.primetest
            && std::io::stderr().is_terminal(),
        last: 0,
        flags: 0,
    };
    let from_stdin = cli.numbers.is_empty() || cli.numbers.iter().any(|n| n == "-");
    for arg in &cli.numbers {
        if arg != "-" && !run.number(arg) {
            return run.exit_code();
        }
    }
    if from_stdin {
        for input in read_stdin() {
            match input {
                Ok(line) => {
                    if !run.number(&line) {
                        break;
                    }
                }
                Err(e) => {
                    eprintln!("ecm-rs: cannot read the standard input: {e}");
                    run.flags |= status::ERROR;
                    break;
                }
            }
        }
    }
    run.exit_code()
}

/// The [`Factorizer`] of the options.
fn factorizer(cli: &Cli) -> Result<Factorizer, String> {
    let param = match (cli.param, &cli.sigma) {
        (Some(p), Some((Some(q), _))) if p != *q => {
            return Err(format!("--param {p} and --sigma {q}:... disagree"));
        }
        (Some(p), _) | (None, &Some((Some(p), _))) => Some(p),
        _ => None,
    };
    let mut f = Factorizer::new();
    if cli.pm1 {
        f = f.algorithm(Algorithm::Pm1);
    }
    if let Some(p) = param {
        f = f.param(Param::try_from(p).map_err(|e| e.to_string())?);
    }
    if let Some(seed) = cli.seed {
        f = f.seed(seed);
    }
    if let Some(b1) = cli.b1 {
        f = f.b1(b1);
    }
    if let Some(b2) = cli.b2 {
        f = f.b2(b2);
    }
    match (cli.curves, &cli.sigma) {
        (Some(curves), _) => f = f.curves(curves),
        // One curve, as GMP-ECM.
        (None, Some(_)) => f = f.curves(1),
        _ => {}
    }
    if let Some((_, sigma)) = &cli.sigma {
        f = f.sigma(sigma.clone());
    }
    if let Some(mb) = cli.maxmem {
        f = f.max_memory(mb.saturating_mul(1 << 20));
    }
    if let Some(timeout) = cli.timeout {
        f = f.timeout(timeout);
    }
    // The options are checked before any work: factoring 1 is immediate.
    match f.clone().factor_partial(&Integer::from(1)).error {
        Some(Error::InvalidOption(msg)) => Err(format!("invalid options: {msg}")),
        Some(Error::BoundsTooSmall) => {
            Err("invalid options: B1 must be at least 6, B2 at least 4".into())
        }
        Some(e) => Err(format!("invalid options: {e}")),
        None => Ok(f),
    }
}

/// Reads the numbers of the standard input: blank lines skipped, `//` comments, `\` line
/// continuation. Invalid UTF-8 is replaced (the expression is then rejected, not the input).
fn read_stdin() -> impl Iterator<Item = std::io::Result<String>> {
    let mut stdin = std::io::stdin().lock();
    let mut next_line = move || -> Option<std::io::Result<String>> {
        let mut bytes = Vec::new();
        match stdin.read_until(b'\n', &mut bytes) {
            Ok(0) => None,
            Ok(_) => {
                let line = String::from_utf8_lossy(&bytes);
                let line = line.strip_suffix('\n').unwrap_or(&line);
                Some(Ok(line.strip_suffix('\r').unwrap_or(line).to_string()))
            }
            Err(e) => Some(Err(e)),
        }
    };
    std::iter::from_fn(move || {
        let mut number = String::new();
        loop {
            let line = match next_line() {
                None if number.trim().is_empty() => return None,
                None => return Some(Ok(number)),
                Some(Err(e)) => return Some(Err(e)),
                Some(Ok(line)) => line,
            };
            let line = line.split_once("//").map_or(line.as_str(), |(l, _)| l);
            let line = line.trim_end();
            match line.strip_suffix('\\') {
                Some(part) => number.push_str(part),
                None => {
                    number.push_str(line);
                    if !number.trim().is_empty() {
                        return Some(Ok(number));
                    }
                    number.clear();
                }
            }
        }
    })
}

/// `input` for a message: its start if it is long.
fn shown(input: &str) -> Cow<'_, str> {
    const MAX: usize = 60;
    match input.char_indices().nth(MAX) {
        Some((end, _)) => format!(
            "{}... ({} characters)",
            &input[..end],
            input.chars().count()
        )
        .into(),
        None => input.into(),
    }
}

/// The factorization of one number, for the output.
struct Outcome {
    primes: Vec<(Integer, usize)>,
    unfactored: Vec<(Integer, usize)>,
    error: Option<Error>,
}

struct Run<'a> {
    cli: &'a Cli,
    factorizer: Factorizer,
    interrupted: Arc<AtomicBool>,
    /// Whether a factorization is in progress (Ctrl-C interrupts it, instead of exiting).
    busy: Arc<AtomicBool>,
    /// End of `--total-timeout`.
    deadline: Option<Instant>,
    progress: bool,
    /// Status bits of the last number.
    last: u8,
    /// Status bits of any number ([`status::ERROR`], [`status::TIMEOUT`]).
    flags: u8,
}

impl Run<'_> {
    fn exit_code(&self) -> ExitCode {
        if self.interrupted.load(Ordering::SeqCst) {
            return ExitCode::from(status::INTERRUPTED);
        }
        ExitCode::from(self.last | self.flags)
    }

    /// Writes a line on stdout (flushed: a pipe sees each result at once).
    fn print(&mut self, line: &str) -> bool {
        let mut out = std::io::stdout().lock();
        if writeln!(out, "{line}").and_then(|()| out.flush()).is_err() {
            self.flags |= status::ERROR;
            return false;
        }
        true
    }

    /// Processes one input: `false` to stop (Ctrl-C, or stdout closed).
    fn number(&mut self, input: &str) -> bool {
        let input = input.trim();
        let n = match expr::eval(input) {
            Ok(n) if n > 0 => n,
            Ok(_) => return self.invalid(input, "the number must be positive".into()),
            Err(e) => return self.invalid(input, format!("invalid expression: {e}")),
        };
        if self.cli.primetest {
            let prime = is_prime(&n);
            self.last = if prime { status::PRIME_COFACTOR } else { 0 };
            return if self.cli.json {
                let value =
                    json!({"input": input, "n": n.to_string(), "prime": prime, "error": null});
                self.print(&value.to_string())
            } else {
                let kind = if prime { "prime" } else { "composite" };
                self.print(&format!("{n}: {kind}"))
            };
        }

        let mut ui = Ui::new(
            if self.cli.quiet { 0 } else { self.cli.verbose },
            self.progress,
            self.cli.b1.is_none() || self.cli.pm1,
        );
        if self.cli.verbose > 0 && !self.cli.quiet {
            ui.line(&format!("Input number is {input} ({} digits)", digits(&n)));
        }
        let start = Instant::now();
        let mut factorizer = self.factorizer.clone();
        if let Some(deadline) = self.deadline {
            let left = deadline.saturating_duration_since(start);
            factorizer = factorizer.timeout(self.cli.timeout.map_or(left, |t| t.min(left)));
        }
        self.busy.store(true, Ordering::SeqCst);
        let outcome = {
            let mut f = factorizer.on_event(|event| {
                ui.handle(event);
                ControlFlow::Continue(())
            });
            if self.cli.one {
                one(&mut f, &n)
            } else {
                let result = f.factor_partial(&n);
                Outcome {
                    primes: result.primes.into_iter().collect(),
                    unfactored: result.unfactored,
                    error: result.error,
                }
            }
        };
        self.busy.store(false, Ordering::SeqCst);
        let time = start.elapsed();
        ui.finish();
        let outcome = normalize(outcome);
        // Whether a proper factor was found: the parts are not just `n`.
        let whole = (n.clone(), 1);
        let split = n > 1
            && !matches!(
                (&outcome.primes[..], &outcome.unfactored[..]),
                ([p], []) | ([], [p]) if *p == whole
            );

        let complete = outcome.unfactored.is_empty();
        self.last = bit(split, status::FACTOR) | bit(complete, status::PRIME_COFACTOR);
        if self.cli.one {
            self.last = one_status(&outcome, split);
        } else if !outcome.primes.is_empty() && split {
            self.last |= status::PRIME_FACTOR;
        }
        let error = outcome.error.as_ref().map(|e| self.error_name(e));
        match error {
            Some("timeout") => {
                self.flags |= status::TIMEOUT;
                if !self.cli.quiet && !self.cli.json {
                    eprintln!(
                        "ecm-rs: timeout after {:.1}s on {}",
                        time.as_secs_f64(),
                        shown(input)
                    );
                }
            }
            Some("interrupted") if !self.cli.quiet && !self.cli.json => {
                eprintln!("ecm-rs: interrupted");
            }
            _ => {}
        }

        let printed = if self.cli.json {
            let list = |parts: &[(Integer, usize)], key: &str| -> Value {
                parts
                    .iter()
                    .map(|(p, e)| json!({key: p.to_string(), "exponent": e}))
                    .collect()
            };
            let value = json!({
                "input": input,
                "n": n.to_string(),
                "factors": list(&outcome.primes, "p"),
                "unfactored": list(&outcome.unfactored, "n"),
                "complete": complete,
                "error": error,
                "time": time.as_secs_f64(),
            });
            self.print(&value.to_string())
        } else {
            let mut parts: Vec<(&Integer, usize, bool)> = outcome
                .primes
                .iter()
                .map(|(p, e)| (p, *e, true))
                .chain(outcome.unfactored.iter().map(|(c, e)| (c, *e, false)))
                .collect();
            parts.sort();
            let mut line = format!("{n} =");
            for (i, (p, e, prime)) in parts.into_iter().enumerate() {
                line.push_str(if i == 0 { " " } else { " * " });
                line.push_str(&product([(p, e)]));
                if !prime {
                    line.push_str(" (composite)");
                }
            }
            if n == 1 {
                line.push_str(" 1");
            }
            self.print(&line)
        };
        printed && !self.interrupted.load(Ordering::SeqCst)
    }

    /// Reports an invalid input: `false` if the output is closed.
    fn invalid(&mut self, input: &str, msg: String) -> bool {
        self.flags |= status::ERROR;
        if self.cli.json {
            let value = if self.cli.primetest {
                json!({"input": input, "n": null, "prime": null, "error": "invalid", "message": msg})
            } else {
                json!({
                    "input": input,
                    "n": null,
                    "factors": [],
                    "unfactored": [],
                    "complete": false,
                    "error": "invalid",
                    "message": msg,
                    "time": 0.0,
                })
            };
            return self.print(&value.to_string());
        }
        eprintln!("ecm-rs: {}: {msg}", shown(input));
        true
    }

    fn error_name(&self, error: &Error) -> &'static str {
        match error {
            Error::Interrupted if self.interrupted.load(Ordering::SeqCst) => "interrupted",
            Error::Interrupted => "timeout",
            Error::ECMFailed => "failed",
            _ => "error",
        }
    }
}

/// `--one`: a factor of `n` and its cofactor.
fn one<H: ecm::EventHandler>(f: &mut Factorizer<H>, n: &Integer) -> Outcome {
    let whole = |error| Outcome {
        primes: Vec::new(),
        unfactored: vec![(n.clone(), 1)],
        error,
    };
    if *n == 1 {
        return Outcome {
            unfactored: Vec::new(),
            ..whole(None)
        };
    }
    match f.find_factor(n) {
        Ok(factor) => {
            let cofactor = Integer::from(n / &factor);
            let (mut primes, mut unfactored) = (Vec::new(), Vec::new());
            for part in [factor, cofactor] {
                if is_prime(&part) {
                    primes.push((part, 1));
                } else {
                    unfactored.push((part, 1));
                }
            }
            Outcome {
                primes,
                unfactored,
                error: None,
            }
        }
        Err(Error::NumberIsPrime) => Outcome {
            primes: vec![(n.clone(), 1)],
            unfactored: Vec::new(),
            error: None,
        },
        Err(error) => whole(Some(error)),
    }
}

/// Merges equal parts, moves the prime parts not tested yet (after an interruption) to the
/// primes, sorts.
fn normalize(mut outcome: Outcome) -> Outcome {
    let mut primes: HashMap<Integer, usize> = HashMap::new();
    let mut unfactored: HashMap<Integer, usize> = HashMap::new();
    for (p, e) in outcome.primes.drain(..) {
        *primes.entry(p).or_default() += e;
    }
    for (c, e) in outcome.unfactored.drain(..) {
        if is_prime(&c) {
            *primes.entry(c).or_default() += e;
        } else {
            *unfactored.entry(c).or_default() += e;
        }
    }
    outcome.primes = primes.into_iter().collect();
    outcome.primes.sort();
    outcome.unfactored = unfactored.into_iter().collect();
    outcome.unfactored.sort();
    outcome
}

/// The status bits of `--one`: GMP-ECM's.
fn one_status(outcome: &Outcome, split: bool) -> u8 {
    if !split {
        return if outcome.unfactored.is_empty() {
            status::PRIME_COFACTOR
        } else {
            0
        };
    }
    // The two parts: the smaller is reported as the factor (as GMP-ECM usually finds it).
    let mut parts: Vec<(&Integer, bool)> = outcome
        .primes
        .iter()
        .map(|(p, _)| (p, true))
        .chain(outcome.unfactored.iter().map(|(c, _)| (c, false)))
        .collect();
    parts.sort();
    let factor_prime = parts.first().is_some_and(|&(_, prime)| prime);
    let cofactor_prime = match parts.as_slice() {
        // p^2: a prime factor, a prime cofactor.
        [_] => factor_prime,
        [_, (_, prime)] => *prime,
        _ => false,
    };
    status::FACTOR
        | bit(factor_prime, status::PRIME_FACTOR)
        | bit(cofactor_prime, status::PRIME_COFACTOR)
}

/// `bit` if `cond`, else 0.
fn bit(cond: bool, bit: u8) -> u8 {
    if cond { bit } else { 0 }
}

/// `--printconfig`.
fn print_config() {
    use gmp_mpfr_sys::gmp;
    println!("ecm-rs {}", env!("CARGO_PKG_VERSION"));
    println!(
        "GMP {}.{}.{} ({}-bit limbs)",
        gmp::VERSION,
        gmp::VERSION_MINOR,
        gmp::VERSION_PATCHLEVEL,
        gmp::LIMB_BITS
    );
    println!(
        "Target: {}-{}",
        std::env::consts::ARCH,
        std::env::consts::OS
    );
    #[cfg(target_arch = "x86_64")]
    println!(
        "BMI2/ADX: {}",
        if config::bmi2_adx_detected() {
            "detected (the stage 1 ladder and the stage 2 pairs use mulx/adcx/adox)"
        } else {
            "not detected (generic code)"
        }
    );
    #[cfg(not(target_arch = "x86_64"))]
    println!("BMI2/ADX: not applicable (not x86_64)");
    println!("Modular arithmetic by size of n (odd n):");
    let (own, max) = (
        if config::MPN_ENABLED {
            config::GMP_LIMBS - 1
        } else {
            config::MAX_LIMBS
        },
        config::MAX_LIMBS,
    );
    let row = |range: String, what: &str| println!("  {range:<20}{what}");
    row(
        format!("up to {} bits", 64 * own),
        "Montgomery, fixed-size limb arrays (own CIOS code)",
    );
    if own < max {
        row(
            format!("{} to {} bits", 64 * own + 1, 64 * max),
            "Montgomery, fixed-size limb arrays (GMP mpn product + REDC)",
        );
    }
    row(
        format!("above {} bits", 64 * max),
        "GMP integers (mpz), plain reduction",
    );
    println!(
        "Default parametrization: {} (GMP-ECM -param)",
        Param::default()
    );
    println!(
        "Stage 2 memory (default): {} MiB",
        config::MAX_POLY_MEMORY >> 20
    );
}
