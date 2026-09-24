//! The `ecm-rs` command line tool (`--features cli`).

use serde_json::Value;
use std::{
    io::Write,
    process::{Command, Output, Stdio},
    time::{Duration, Instant},
};

/// A number with a 25-digit factor: out of reach of a short timeout.
const P25: &str = "158200595350375547764951423818039293423703269411061670349759";

fn ecm_rs(args: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_ecm-rs"));
    command.args(args);
    command
}

/// Runs `ecm-rs args` with `stdin`: (exit code, stdout, stderr).
fn run_with(args: &[&str], stdin: &str) -> (i32, String, String) {
    let mut child = ecm_rs(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(stdin.as_bytes())
        .unwrap();
    let Output {
        status,
        stdout,
        stderr,
    } = child.wait_with_output().unwrap();
    (
        status.code().unwrap(),
        String::from_utf8(stdout).unwrap(),
        String::from_utf8(stderr).unwrap(),
    )
}

fn run(args: &[&str]) -> (i32, String, String) {
    run_with(args, "")
}

#[test]
fn factors_completely() {
    let (code, out, err) = run(&["4516511326451341281684513", "17", "1"]);
    assert_eq!(
        out,
        "4516511326451341281684513 = 3^2 * 39869 * 131743543 * 95542348571\n17 = 17\n1 = 1\n"
    );
    assert_eq!(err, "");
    // The last number is 1: nothing to factor.
    assert_eq!(code, 8);
    let (code, out, _) = run(&["7060005655815754299976961394452809"]);
    assert_eq!(
        out,
        "7060005655815754299976961394452809 = 6988699669998001 * 1010203040506070809\n"
    );
    assert_eq!(code, 14);
}

#[test]
fn expressions() {
    let (code, out, _) = run(&["2^67-1", "10!+1", "(10^20-1)/9", "30#+1", "2^64", "17#5"]);
    assert_eq!(
        out,
        "147573952589676412927 = 193707721 * 761838257287\n\
         3628801 = 11 * 329891\n\
         11111111111111111111 = 11 * 41 * 101 * 271 * 3541 * 9091 * 27961\n\
         6469693231 = 331 * 571 * 34231\n\
         18446744073709551616 = 2^64\n\
         85085 = 5 * 7 * 11 * 13 * 17\n"
    );
    assert_eq!(code, 14);
}

#[test]
fn stdin_lines() {
    let input = "398883434337287\n\n  // a comment\n2^67\\\n-1 // continued\n   \n1000003\n";
    let (code, out, _) = run_with(&[], input);
    assert_eq!(
        out,
        "398883434337287 = 4009823 * 99476569\n\
         147573952589676412927 = 193707721 * 761838257287\n\
         1000003 = 1000003\n"
    );
    assert_eq!(code, 8);
    // "-" reads stdin after the arguments before it.
    let (code, out, _) = run_with(&["15", "-"], "21\n");
    assert_eq!(out, "15 = 3 * 5\n21 = 3 * 7\n");
    assert_eq!(code, 14);
}

#[test]
fn invalid_inputs() {
    let (code, out, err) = run(&["12a", "0", "(1+2", "7/2", "10^10^10", "100000!", "6"]);
    assert_eq!(out, "6 = 2 * 3\n");
    for msg in [
        "12a: invalid expression: unexpected 'a' (at character 3)",
        "0: the number must be positive",
        "(1+2: invalid expression: unbalanced '(' (at character 1)",
        "7/2: invalid expression: inexact division (at character 2)",
        "10^10^10: invalid expression: value too large",
        "100000!: invalid expression: argument of '!' must be between 0 and 50000",
    ] {
        assert!(err.contains(msg), "{msg} not in {err}");
    }
    // Error bit, and the result of the last number.
    assert_eq!(code, 1 | 14);

    let (code, out, _) = run(&["--json", "--", "-5"]);
    let value: Value = serde_json::from_str(out.trim()).unwrap();
    assert_eq!(value["error"], "the number must be positive");
    assert_eq!(value["n"], Value::Null);
    assert_eq!(code, 1);
}

#[test]
fn usage_errors() {
    for args in [
        &["--b2", "1000", "15"][..],
        &["-c", "10", "15"],
        &["--sigma", "5", "15"],
        &["--param", "3", "15"],
        &["--b1", "abc", "15"],
        &["--timeout", "-1", "15"],
        &["--unknown"],
        &["--b1", "11000", "--sigma", "1:5", "--param", "2", "15"],
        &["--b1", "11000", "--sigma", "1:1", "15"],
        &["--b1", "4", "15"],
        &["--pm1", "--b1", "1000", "--sigma", "7", "15"],
        &["--primetest", "--one", "15"],
        &["-q", "-v", "15"],
    ] {
        let (code, out, err) = run(args);
        assert_eq!(code, 64, "{args:?}: {err}");
        assert_eq!(out, "", "{args:?}");
        assert!(!err.is_empty(), "{args:?}");
    }
    let (code, out, _) = run(&["--help"]);
    assert_eq!(code, 0);
    assert!(out.contains("GMP-ECM equivalents"));
    assert!(out.contains("Exit status"));
    let (code, out, _) = run(&["--version"]);
    assert_eq!(code, 0);
    assert_eq!(out, format!("ecm-rs {}\n", env!("CARGO_PKG_VERSION")));
}

#[test]
fn json() {
    let (code, out, err) = run(&["--json", "-v", "2^67-1", "3^5*1000003"]);
    assert_eq!(code, 14);
    assert!(err.contains("Factor found"), "-v still on stderr");
    let lines: Vec<Value> = out
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(lines.len(), 2);
    let v = &lines[0];
    assert_eq!(v["input"], "2^67-1");
    assert_eq!(v["n"], "147573952589676412927");
    assert_eq!(
        v["factors"],
        serde_json::json!([
            {"p": "193707721", "exponent": 1},
            {"p": "761838257287", "exponent": 1},
        ])
    );
    assert_eq!(v["unfactored"], serde_json::json!([]));
    assert_eq!(v["complete"], true);
    assert_eq!(v["error"], Value::Null);
    assert!(v["time"].as_f64().unwrap() >= 0.0);
    assert_eq!(
        lines[1]["factors"],
        serde_json::json!([{"p": "3", "exponent": 5}, {"p": "1000003", "exponent": 1}])
    );
}

#[test]
fn one_factor() {
    let (code, out, _) = run(&["--one", "2^67-1"]);
    assert_eq!(out, "147573952589676412927 = 193707721 * 761838257287\n");
    assert_eq!(code, 14);
    // Trial division finds 2: the cofactor is composite.
    let (code, out, _) = run(&["--one", "1000004"]);
    assert_eq!(out, "1000004 = 2 * 500002 (composite)\n");
    assert_eq!(code, 6);
    let (code, out, _) = run(&["--one", "1000003"]);
    assert_eq!(out, "1000003 = 1000003\n");
    assert_eq!(code, 8);
    let (code, out, _) = run(&["--one", "--json", "1000004"]);
    let v: Value = serde_json::from_str(out.trim()).unwrap();
    assert_eq!(v["complete"], false);
    assert_eq!(
        v["unfactored"],
        serde_json::json!([{"n": "500002", "exponent": 1}])
    );
    assert_eq!(code, 6);
}

#[test]
fn pm1_only() {
    // 193707721 - 1 = 2^3 * 3^3 * 5 * 11 * 11 * 107 * ... smooth enough for B1 = 6000.
    let (code, out, err) = run(&["--pm1", "--b1", "6000", "-v", "2^67-1"]);
    assert_eq!(out, "147573952589676412927 = 193707721 * 761838257287\n");
    assert!(err.contains("P-1 on C21: B1=6000"), "{err}");
    assert!(err.contains("Factor found by P-1 stage"), "{err}");
    assert!(!err.contains("Curve"), "{err}");
    assert_eq!(code, 14);
    // Too small bounds: nothing found.
    let (code, out, _) = run(&["--pm1", "--b1", "10", "--b2", "10", "2^67-1"]);
    assert_eq!(
        out,
        "147573952589676412927 = 147573952589676412927 (composite)\n"
    );
    assert_eq!(code, 0);
}

#[test]
fn sigma_reproduces_a_curve() {
    // GMP-ECM 7.0.7: `echo N | ecm -sigma 1:1176292814 11000` finds the 20-digit factor in
    // step 2.
    let n = "15658598057181786459081452046251445462002559800474409088109";
    let expected =
        format!("{n} = 13507140964289979319 * 1159282937712710938601499347662537052411\n");
    let (code, out, err) = run(&["--b1", "11000", "--sigma", "1:1176292814", "-v", n]);
    assert_eq!(out, expected);
    assert!(
        err.contains("Curve 1/1: sigma=1:1176292814"),
        "one curve by default: {err}"
    );
    assert!(err.contains("Factor found by ECM stage 2"), "{err}");
    assert_eq!(code, 14);
    // Same with --param 1 and a bare sigma.
    let (_, out, _) = run(&["--b1", "11000", "--param", "1", "--sigma", "1176292814", n]);
    assert_eq!(out, expected);
    // The previous sigma misses it.
    let (code, out, _) = run(&["--b1", "11000", "--sigma", "1:1176292813", n]);
    assert_eq!(out, format!("{n} = {n} (composite)\n"));
    assert_eq!(code, 0);
}

/// The `sigma=...` of the `-v` curve lines.
fn sigmas(err: &str) -> Vec<&str> {
    err.lines()
        .filter(|line| line.starts_with("Curve "))
        .map(|line| {
            line.split(", ")
                .next()
                .unwrap()
                .split(' ')
                .next_back()
                .unwrap()
        })
        .collect()
}

#[test]
fn seed_determinism() {
    let n = "15658598057181786459081452046251445462002559800474409088109";
    let args = |seed| ["--b1", "2000", "-c", "3", "-v", "--seed", seed, n];
    let (_, out1, err1) = run(&args("7"));
    let (_, out2, err2) = run(&args("7"));
    let (_, _, err3) = run(&args("8"));
    assert_eq!(out1, out2);
    assert_eq!(sigmas(&err1).len(), 3);
    assert_eq!(sigmas(&err1), sigmas(&err2));
    assert_ne!(sigmas(&err1), sigmas(&err3));
    // The default seed is fixed too.
    let (_, _, a) = run(&["--b1", "2000", "-c", "2", "-v", n]);
    let (_, _, b) = run(&["--b1", "2000", "-c", "2", "-v", n]);
    assert_eq!(sigmas(&a), sigmas(&b));
}

#[test]
fn timeout_prints_partial_results() {
    let n = format!("3*1000003*{P25}");
    let start = Instant::now();
    let (code, out, err) = run(&["--timeout", "0.3", &n, "15"]);
    assert!(start.elapsed() < Duration::from_secs(5));
    let product = P25.parse::<rug::Integer>().unwrap() * 3_000_009u32;
    assert!(
        out.starts_with(&format!("{product} = 3 * 1000003 * {P25} (composite)\n")),
        "{out}"
    );
    // The next number is still factored.
    assert!(out.ends_with("15 = 3 * 5\n"), "{out}");
    assert!(err.contains("timeout"), "{err}");
    assert_eq!(code, 16 | 14);

    let (code, out, _) = run(&["--timeout", "0.2", "--json", P25]);
    let v: Value = serde_json::from_str(out.trim()).unwrap();
    assert_eq!(v["error"], "timeout");
    assert_eq!(v["complete"], false);
    assert_eq!(v["factors"], serde_json::json!([]));
    assert_eq!(code, 16);
}

#[cfg(unix)]
#[test]
fn ctrl_c_prints_partial_results() {
    let child = ecm_rs(&[&format!("5*{P25}"), "15"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    std::thread::sleep(Duration::from_millis(500));
    let killed = Command::new("kill")
        .args(["-INT", &child.id().to_string()])
        .status()
        .unwrap();
    assert!(killed.success());
    let output = child.wait_with_output().unwrap();
    let out = String::from_utf8(output.stdout).unwrap();
    assert!(
        out.ends_with(&format!(" = 5 * {P25} (composite)\n")),
        "{out}"
    );
    // The next numbers are not processed.
    assert_eq!(out.lines().count(), 1);
    assert_eq!(output.status.code(), Some(130));
}

#[test]
fn primetest() {
    let (code, out, _) = run(&["--primetest", "2^61-1", "2^67-1"]);
    assert_eq!(
        out,
        "2305843009213693951: prime\n147573952589676412927: composite\n"
    );
    assert_eq!(code, 0);
    let (code, out, _) = run(&["--primetest", "--json", "2^61-1"]);
    let v: Value = serde_json::from_str(out.trim()).unwrap();
    assert_eq!(v["prime"], true);
    assert_eq!(code, 8);
}

#[test]
fn printconfig() {
    let (code, out, _) = run(&["--printconfig"]);
    assert_eq!(code, 0);
    assert!(out.starts_with(&format!("ecm-rs {}\n", env!("CARGO_PKG_VERSION"))));
    assert!(out.contains("GMP "));
    assert!(out.contains("BMI2/ADX: "));
    assert!(out.contains("Montgomery"));
    let (code, _, _) = run(&["--printconfig", "15"]);
    assert_eq!(code, 64);
}

#[test]
fn quiet_and_bounds() {
    let (code, out, err) = run(&["-q", "--b1", "11e3", "-c", "50", "2^67-1"]);
    assert_eq!(out, "147573952589676412927 = 193707721 * 761838257287\n");
    assert_eq!(err, "");
    assert_eq!(code, 14);
    // Odd bounds are rounded up to even (the same primes).
    let (_, out, err) = run(&["-v", "--b1", "2001", "--b2", "100001", "-c", "1", "2^67-1"]);
    assert!(err.contains("Using B1=2002, B2="), "{err}");
    assert!(!out.is_empty());
}
