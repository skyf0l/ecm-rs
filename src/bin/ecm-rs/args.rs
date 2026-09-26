//! The command line: options, their parsing and the help texts.

use clap::Parser;
use rug::Integer;
use std::time::Duration;

/// Exit status bits, as GMP-ECM's.
pub mod status {
    /// An input was invalid, or the output could not be written.
    pub const ERROR: u8 = 1;
    /// A proper factor was found.
    pub const FACTOR: u8 = 2;
    /// A prime factor was found (with `--one`: the factor found is prime).
    pub const PRIME_FACTOR: u8 = 4;
    /// The cofactor is prime or 1: the factorization is complete (with `--one`: the cofactor
    /// of the factor found is prime; with `--primetest`: the number is prime).
    pub const PRIME_COFACTOR: u8 = 8;
    /// `--timeout` or `--total-timeout` interrupted a factorization.
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
\"(composite)\" when the factorization is incomplete (--one, --timeout, Ctrl-C, --curves, \
--pm1 or --pp1 exhausted). Diagnostics, -v lines and the progress bar go to stderr.";

const AFTER_LONG_HELP: &str = "\
GMP-ECM equivalents:
  ecm B1                       ecm-rs --b1 B1 -c 1     (GMP-ECM runs one curve by default)
  ecm B1 B2                    ecm-rs --b1 B1 --b2 B2 -c 1
  ecm -c N B1                  ecm-rs --b1 B1 -c N     (N curves per composite part)
  ecm -sigma 1:S B1            ecm-rs --b1 B1 --sigma 1:S   (one curve, as GMP-ECM)
  ecm -param 0 B1              ecm-rs --b1 B1 -c 1 --param 0
  ecm -one ...                 ecm-rs --one ...
  ecm -pm1 B1 B2               ecm-rs --pm1 --b1 B1 --b2 B2
  ecm -pp1 -x0 2/7 B1 B2       ecm-rs --pp1 --b1 B1 --b2 B2   (x0 = 2/7 by default; GMP-ECM:
                                 random, which finds a smooth p+1 half of the time)
  ecm -pm1/-pp1 -x0 X ...      ecm-rs --pm1/--pp1 --x0 X ...
  ecm -maxmem MB               ecm-rs --maxmem MB
  ecm -base2 K / -nobase2      ecm-rs --base2 K / --nobase2   (by default, as GMP-ECM: the
                                 numbers dividing 2^k+-1 compute modulo it, when faster)
  ecm -primetest               ecm-rs --primetest
  ecm -printconfig             ecm-rs --printconfig
  ecm -q / -v                  ecm-rs -q / -v
  (no equivalent)              ecm-rs --b1 B1    (curves until a factor is found)
  (no equivalent)              ecm-rs N    (complete factorization, bounds by factor size)

Exit status (bits, as GMP-ECM's, for the last number; 1 and 16 for any number):
  0    no factor found (curves exhausted, or --primetest: composite)
  1    error: an invalid number (the others are still processed)
  2    a composite factor found, the cofactor is composite (--one)
  6    a prime factor found, the cofactor is composite (incomplete factorization)
  8    the number is prime (or 1): nothing to factor
  10   a composite factor found, the cofactor is prime (--one)
  14   factored completely (--one: a prime factor and a prime cofactor)
  +16  --timeout or --total-timeout interrupted a factorization (its partial result is printed)
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
pub struct Cli {
    /// Numbers (or expressions) to factor; with none, or "-", read from stdin, one per line.
    #[arg(value_name = "N")]
    pub numbers: Vec<String>,

    /// Fixed stage 1 bound (e.g. 11000, 11e3): curves with these bounds instead of the bounds
    /// by factor size.
    #[arg(long, value_name = "B1", value_parser = parse_bound)]
    pub b1: Option<usize>,

    /// Stage 2 bound [default: GMP-ECM's default for B1].
    #[arg(long, value_name = "B2", value_parser = parse_bound, requires = "b1")]
    pub b2: Option<usize>,

    /// Largest number of curves per composite part [default: unlimited; 1 with --sigma].
    #[arg(short, long, value_name = "N", requires = "b1", conflicts_with_all = ["pm1", "pp1"])]
    pub curves: Option<usize>,

    /// Parameter of the first curve, "P:S" or "S", with P the parametrization (the next curves
    /// take S+1, S+2...). Reproduces a curve of GMP-ECM ("-sigma P:S") or of -v. S is at least
    /// 6 with P = 0, in [2, 2^32) with P = 1, in [2, 2^64) with P = 2.
    #[arg(
        long,
        value_name = "[P:]S",
        value_parser = parse_sigma,
        requires = "b1",
        conflicts_with_all = ["pm1", "pp1"]
    )]
    pub sigma: Option<(Option<u8>, Integer)>,

    /// Parametrization of the curves, as GMP-ECM's -param: 0 (Suyama), 1, 2 [default: 2].
    #[arg(
        long,
        value_name = "P",
        value_parser = clap::value_parser!(u8).range(0..=2),
        conflicts_with_all = ["pm1", "pp1"]
    )]
    pub param: Option<u8>,

    /// Stops at the first factor found (maybe composite): prints it and its cofactor.
    #[arg(long)]
    pub one: bool,

    /// Runs only Pollard's P-1 method (with --b1: once per composite part).
    #[arg(long, conflicts_with = "pp1")]
    pub pm1: bool,

    /// Runs only Williams' P+1 method (with --b1: once per composite part): finds the factors
    /// p with a smooth p+1 (or p-1, depending on p and the seed).
    #[arg(long)]
    pub pp1: bool,

    /// Starting value of P-1 or P+1, an integer or a fraction "N/D" (as GMP-ECM's -x0)
    /// [default: 3 for P-1, 2/7 for P+1].
    #[arg(long, value_name = "X", value_parser = parse_x0, allow_hyphen_values = true)]
    pub x0: Option<(Integer, Integer)>,

    /// Seed of the random curves [default: fixed, the results are reproducible].
    #[arg(long, value_name = "SEED")]
    pub seed: Option<u64>,

    /// Threads running the curves (and P-1 alongside them) [default, or 0: the available
    /// parallelism, at most 1024]. The results (factors, curves, -v lines but for the times) do
    /// not depend on it.
    #[arg(short, long, value_name = "N", value_parser = clap::value_parser!(u32).range(..=1024))]
    pub threads: Option<u32>,

    /// Gives up on each number after SECS seconds, printing what was found.
    #[arg(long, value_name = "SECS", value_parser = parse_timeout)]
    pub timeout: Option<Duration>,

    /// Gives up after SECS seconds in total: the number being factored and the next ones are
    /// printed with what was found (at least the small factors of trial division).
    #[arg(long, value_name = "SECS", value_parser = parse_timeout)]
    pub total_timeout: Option<Duration>,

    /// Memory for stage 2, in MiB, per thread [default: 256].
    #[arg(long, value_name = "MB")]
    pub maxmem: Option<usize>,

    /// Computes modulo 2^K+1 (K > 0) or 2^-K-1 (K < 0), with a reduction by shifts and
    /// additions, as GMP-ECM's -base2: every composite part must divide it. By default, the
    /// numbers that divide 2^k+-1 with k at most 1.4 times their size do when it is faster.
    #[arg(
        long,
        value_name = "K",
        allow_hyphen_values = true,
        value_parser = parse_base2,
        conflicts_with = "nobase2"
    )]
    pub base2: Option<i64>,

    /// Never computes modulo 2^k+-1 (GMP-ECM's -nobase2).
    #[arg(long)]
    pub nobase2: bool,

    /// Only tests whether each number is prime: prints "N: prime" or "N: composite".
    #[arg(long, conflicts_with_all = ["one", "pm1", "pp1", "b1"])]
    pub primetest: bool,

    /// Prints the configuration (versions, arithmetic by size, CPU extensions) and exits.
    #[arg(long, exclusive = true)]
    pub printconfig: bool,

    /// Only the results: no diagnostics, no progress bar.
    #[arg(short, long, conflicts_with = "verbose")]
    pub quiet: bool,

    /// Prints the steps (as GMP-ECM's -v) on stderr: levels, curves with their sigma and
    /// stage times, P-1 and P+1 runs, factors found. Twice: also the prime factors as they are
    /// found.
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub verbose: u8,

    /// One JSON object per number, on one line: {"input", "n", "factors": [{"p", "exponent"}],
    /// "unfactored": [{"n", "exponent"}], "complete", "error", "time"} (numbers as strings,
    /// time in seconds; "error": null, "timeout", "interrupted", "failed" or "invalid", with
    /// a "message" and "n": null for an invalid input; "prime" instead of the factors with
    /// --primetest).
    #[arg(long)]
    pub json: bool,

    /// No progress bar (shown on stderr when it is a terminal).
    #[arg(long)]
    pub no_progress: bool,
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

/// `N` or `N/D` (integers, `D != 0`).
fn parse_x0(s: &str) -> Result<(Integer, Integer), String> {
    let invalid = || format!("invalid x0 '{s}' (an integer or N/D)");
    let (num, den) = s.split_once('/').unwrap_or((s, "1"));
    let num = num.trim().parse::<Integer>().map_err(|_| invalid())?;
    let den = den.trim().parse::<Integer>().map_err(|_| invalid())?;
    if den == 0 {
        return Err(format!("invalid x0 '{s}': zero denominator"));
    }
    Ok((num, den))
}

/// A non-zero exponent `K` of `--base2`.
fn parse_base2(s: &str) -> Result<i64, String> {
    s.parse::<i64>()
        .ok()
        .filter(|&k| k != 0)
        .ok_or_else(|| format!("invalid exponent '{s}' (K for 2^K+1, -K for 2^K-1)"))
}

fn parse_timeout(s: &str) -> Result<Duration, String> {
    s.parse::<f64>()
        .ok()
        .and_then(|secs| Duration::try_from_secs_f64(secs).ok())
        .filter(|d| !d.is_zero())
        .ok_or_else(|| format!("invalid timeout '{s}' (positive seconds)"))
}
