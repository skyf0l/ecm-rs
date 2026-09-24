//! `ecm-rs`: factors integers with the `ecm` crate, with GMP-ECM-like options.

mod expr;
mod ui;

use clap::{Parser, error::ErrorKind};
use ecm::{Error, Factorizer, Param};
use rug::Integer;
use serde_json::{Value, json};
use std::{
    collections::HashMap,
    io::{BufRead, IsTerminal, Write},
    ops::ControlFlow,
    process::ExitCode,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    time::{Duration, Instant},
};
use ui::{Ui, digits, is_prime, product};

/// Exit status bits, as GMP-ECM's.
mod status {
    /// An input was invalid, or the output could not be written.
    pub const ERROR: u8 = 1;
    /// A proper factor was found.
    pub const FACTOR: u8 = 2;
    /// A prime factor was found (with `--one`: the factor found is prime).
    pub const PRIME_FACTOR: u8 = 4;
    /// The cofactor is prime or 1: the factorization is complete (with `--one`: the cofactor
    /// of the factor found is prime; with `--primetest`: the number is prime).
    pub const PRIME_COFACTOR: u8 = 8;
    /// `--timeout` interrupted a factorization.
    pub const TIMEOUT: u8 = 16;
    /// Invalid command line.
    pub const USAGE: u8 = 64;
    /// Interrupted by Ctrl-C (128 + SIGINT).
    pub const INTERRUPTED: u8 = 130;
}

const LONG_ABOUT: &str = "\
Factors integers with Lenstra's elliptic curve method (ECM), a port of GMP-ECM's techniques.

By default, each number is factored completely: trial division, then P-1 and curves with the \
bounds for factors of 10, 15, 20, ... digits in turn, so that the factors are found from the \
smallest to the largest. With --b1, curves with fixed bounds are run instead, as GMP-ECM does.

The numbers are the arguments, or the lines of the standard input (with no argument, or \"-\"). \
They can be expressions: integers, + - * / (exact division), ^ (right associative), \
parentheses, n! (factorial), n!k (multi-factorial: 15!3 = 15*12*9*6*3), n# (primorial: \
11# = 2*3*5*7*11) and n#k (reduced primorial: 17#5 = 5*7*11*13*17). On the standard input, \
blank lines are skipped, \"//\" starts a comment and a line ending with \"\\\" continues on \
the next one.

Output (stdout), one line per number: \"N = p1^e1 * p2 * ...\", in increasing order, where N is \
the number in decimal. Every factor is a (probable) prime, except the parts marked \
\"(composite)\" when the factorization is incomplete (--one, --timeout, Ctrl-C, --curves or \
--pm1 exhausted). Diagnostics, -v lines and the progress bar go to stderr.";

const AFTER_LONG_HELP: &str = "\
GMP-ECM equivalents:
  ecm B1                       ecm-rs --b1 B1          (curves until a factor is found)
  ecm B1 B2                    ecm-rs --b1 B1 --b2 B2
  ecm -c N B1                  ecm-rs --b1 B1 -c N
  ecm -sigma 1:S B1            ecm-rs --b1 B1 --sigma 1:S   (one curve, as GMP-ECM)
  ecm -param 0 B1              ecm-rs --b1 B1 --param 0
  ecm -one ...                 ecm-rs --one ...
  ecm -pm1 B1 B2               ecm-rs --pm1 --b1 B1 --b2 B2
  ecm -maxmem MB               ecm-rs --maxmem MB
  ecm -primetest               ecm-rs --primetest
  ecm -printconfig             ecm-rs --printconfig
  ecm -q / -v                  ecm-rs -q / -v
  (no equivalent)              ecm-rs N    (complete factorization, bounds by factor size)

Exit status (bits, as GMP-ECM's, for the last number; 1 and 16 for any number):
  0    no factor found (curves exhausted, or --primetest: composite)
  1    error: an invalid number (the others are still processed)
  2    a composite factor found, the cofactor is composite (--one)
  6    a prime factor found, the cofactor is composite (incomplete factorization)
  8    the number is prime (or 1): nothing to factor
  10   a composite factor found, the cofactor is prime (--one)
  14   factored completely (--one: a prime factor and a prime cofactor)
  +16  --timeout interrupted a factorization (its partial result is printed)
  64   invalid command line
  130  interrupted by Ctrl-C (the partial result of the current number is printed)";

/// Factors integers with the elliptic curve method.
#[derive(Debug, Parser)]
#[command(
    name = "ecm-rs",
    version,
    long_about = LONG_ABOUT,
    after_long_help = AFTER_LONG_HELP,
    after_help = "See --help for the input syntax, GMP-ECM equivalents and exit codes."
)]
struct Cli {
    /// Numbers (or expressions) to factor; with none, or "-", read from stdin, one per line.
    #[arg(value_name = "N")]
    numbers: Vec<String>,

    /// Fixed stage 1 bound (e.g. 11000, 11e3): curves with these bounds instead of the bounds
    /// by factor size.
    #[arg(long, value_name = "B1", value_parser = parse_bound)]
    b1: Option<usize>,

    /// Stage 2 bound [default: GMP-ECM's default for B1].
    #[arg(long, value_name = "B2", value_parser = parse_bound, requires = "b1")]
    b2: Option<usize>,

    /// Largest number of curves per composite part [default: unlimited; 1 with --sigma].
    #[arg(short, long, value_name = "N", requires = "b1")]
    curves: Option<usize>,

    /// Parameter of the first curve, "P:S" or "S", with P the parametrization (the next curves
    /// take S+1, S+2...). Reproduces a curve of GMP-ECM ("-sigma P:S") or of -v.
    #[arg(long, value_name = "[P:]S", value_parser = parse_sigma, requires = "b1")]
    sigma: Option<(Option<u8>, Integer)>,

    /// Parametrization of the curves, as GMP-ECM's -param: 0 (Suyama), 1, 2 [default: 2].
    #[arg(long, value_name = "P", value_parser = clap::value_parser!(u8).range(0..=2))]
    param: Option<u8>,

    /// Stops at the first factor found (maybe composite): prints it and its cofactor.
    #[arg(long)]
    one: bool,

    /// Runs only Pollard's P-1 method (with --b1: once per composite part).
    #[arg(long)]
    pm1: bool,

    /// Seed of the random curves [default: fixed, the results are reproducible].
    #[arg(long, value_name = "SEED")]
    seed: Option<u64>,

    /// Gives up on a number after SECS seconds (each number), printing what was found.
    #[arg(long, value_name = "SECS", value_parser = parse_timeout)]
    timeout: Option<Duration>,

    /// Memory for stage 2, in MiB [default: 256].
    #[arg(long, value_name = "MB")]
    maxmem: Option<usize>,

    /// Only tests whether each number is prime: prints "N: prime" or "N: composite".
    #[arg(long, conflicts_with_all = ["one", "pm1", "b1"])]
    primetest: bool,

    /// Prints the configuration (versions, arithmetic by size, CPU extensions) and exits.
    #[arg(long, exclusive = true)]
    printconfig: bool,

    /// Only the results: no diagnostics, no progress bar.
    #[arg(short, long, conflicts_with = "verbose")]
    quiet: bool,

    /// Prints the steps (as GMP-ECM's -v) on stderr: levels, curves with their sigma and
    /// stage times, P-1 runs, factors found. Twice: also the prime factors as they are found.
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    /// One JSON object per number, on one line: {"input", "n", "factors": [{"p", "exponent"}],
    /// "unfactored": [{"n", "exponent"}], "complete", "error", "time"} (numbers as strings,
    /// time in seconds; "prime" instead of the factors with --primetest).
    #[arg(long)]
    json: bool,

    /// No progress bar (shown on stderr when it is a terminal).
    #[arg(long)]
    no_progress: bool,
}

/// A bound: an integer, or in scientific notation (`11e3`, `1.5e6`); odd bounds are rounded
/// up to even (the same primes).
fn parse_bound(s: &str) -> Result<usize, String> {
    let value = if let Ok(v) = s.parse::<usize>() {
        v
    } else {
        let v: f64 = s
            .parse()
            .map_err(|_| format!("invalid bound '{s}' (e.g. 11000 or 11e3)"))?;
        if !(v.is_finite() && v >= 0.0 && v < 2f64.powi(63) && v.fract() == 0.0) {
            return Err(format!("invalid bound '{s}'"));
        }
        v as usize
    };
    Ok(value + (value & 1))
}

/// `[P:]S`.
fn parse_sigma(s: &str) -> Result<(Option<u8>, Integer), String> {
    let (param, sigma) = match s.split_once(':') {
        Some((p, sigma)) => {
            let p = p
                .parse::<u8>()
                .ok()
                .filter(|&p| p <= 2)
                .ok_or_else(|| format!("unsupported parametrization '{p}' (0, 1 or 2)"))?;
            (Some(p), sigma)
        }
        None => (None, s),
    };
    let sigma = sigma
        .parse::<Integer>()
        .map_err(|_| format!("invalid sigma '{sigma}'"))?;
    Ok((param, sigma))
}

fn parse_timeout(s: &str) -> Result<Duration, String> {
    s.parse::<f64>()
        .ok()
        .and_then(|secs| Duration::try_from_secs_f64(secs).ok())
        .filter(|d| !d.is_zero())
        .ok_or_else(|| format!("invalid timeout '{s}' (positive seconds)"))
}

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
    {
        let flag = interrupted.clone();
        // A second Ctrl-C exits at once. No handler (an error): Ctrl-C kills the process.
        let _ = ctrlc::set_handler(move || {
            if flag.swap(true, Ordering::SeqCst) {
                std::process::exit(status::INTERRUPTED.into());
            }
        });
    }
    let factorizer = factorizer.interrupt_flag(interrupted.clone());

    let mut run = Run {
        cli: &cli,
        factorizer,
        interrupted,
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
    let mut f = Factorizer::new().pm1(cli.pm1);
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
/// continuation.
fn read_stdin() -> impl Iterator<Item = std::io::Result<String>> {
    let mut lines = std::io::stdin().lock().lines();
    std::iter::from_fn(move || {
        let mut number = String::new();
        loop {
            let line = match lines.next() {
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
                self.print(&json!({"input": input, "n": n.to_string(), "prime": prime}).to_string())
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
        let outcome = {
            let mut f = self.factorizer.clone().on_event(|event| {
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
                        "ecm-rs: timeout after {:.1}s on {input}",
                        time.as_secs_f64()
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
            return self.print(&json!({"input": input, "n": null, "error": msg}).to_string());
        }
        eprintln!("ecm-rs: {input}: {msg}");
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
    {
        let detected =
            std::is_x86_feature_detected!("bmi2") && std::is_x86_feature_detected!("adx");
        println!(
            "BMI2/ADX: {}",
            if detected {
                "detected (the hot loops of stage 1 and stage 2 use mulx/adcx/adox)"
            } else {
                "not detected (generic code)"
            }
        );
    }
    #[cfg(not(target_arch = "x86_64"))]
    println!("BMI2/ADX: not applicable (not x86_64)");
    // As src/arith.rs: MAX_LIMBS = 16, GMP_LIMBS = 11, and GMP's mpn functions need 64-bit
    // limbs and GMP >= 5.1.
    let mpn = gmp::LIMB_BITS == 64
        && gmp::NAIL_BITS == 0
        && (gmp::VERSION > 5 || (gmp::VERSION == 5 && gmp::VERSION_MINOR >= 1));
    println!("Modular arithmetic by size of n (odd n):");
    if mpn {
        println!("  up to 640 bits      Montgomery, fixed-size limb arrays (own CIOS code)");
        println!(
            "  641 to 1024 bits    Montgomery, fixed-size limb arrays (GMP mpn product + REDC)"
        );
    } else {
        println!("  up to 1024 bits     Montgomery, fixed-size limb arrays (own CIOS code)");
    }
    println!("  above 1024 bits     GMP integers (mpz), plain reduction");
    println!(
        "Default parametrization: {} (GMP-ECM -param)",
        Param::default()
    );
    println!("Stage 2 memory (default): 256 MiB");
}
