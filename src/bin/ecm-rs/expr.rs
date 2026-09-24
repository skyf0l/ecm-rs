//! Input expressions, as GMP-ECM reads them: integers, `+ - * /`, `^` (right associative),
//! parentheses, factorial `n!` (multi-factorial `n!k`) and primorial `n#` (reduced primorial
//! `n#k`, the product of the primes from `k` to `n`).

use rug::{Complete, Integer, ops::Pow};
use std::fmt;

/// Largest size (in bits) of a value, intermediate or final: far beyond what ECM can factor,
/// but it keeps a typo (`10^10^10`) from exhausting the memory.
pub const MAX_BITS: u32 = 1 << 20;

/// Largest argument of `n!` and `n#` (whose values stay below [`MAX_BITS`]).
const MAX_FACTORIAL: u32 = 50_000;
const MAX_PRIMORIAL: u32 = 700_000;

/// Largest nesting of parentheses, signs and exponents (the parser is recursive).
const MAX_DEPTH: usize = 200;

/// Why an expression was rejected, with the position (in characters, from 1) of the problem.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExprError {
    pub pos: usize,
    pub msg: String,
}

impl fmt::Display for ExprError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (at character {})", self.msg, self.pos)
    }
}

/// Evaluates `input`.
pub fn eval(input: &str) -> Result<Integer, ExprError> {
    let mut parser = Parser {
        chars: input.chars().collect(),
        pos: 0,
        depth: 0,
    };
    parser.skip_spaces();
    if parser.peek().is_none() {
        return Err(parser.error("empty expression"));
    }
    let value = parser.expr()?;
    parser.skip_spaces();
    match parser.peek() {
        None => Ok(value),
        Some(')') => Err(parser.error("unbalanced ')'")),
        Some(c) => Err(parser.error(format!("unexpected '{c}'"))),
    }
}

struct Parser {
    chars: Vec<char>,
    pos: usize,
    /// Nesting of [`Parser::unary`], which every recursion goes through.
    depth: usize,
}

impl Parser {
    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn skip_spaces(&mut self) {
        while self.peek().is_some_and(char::is_whitespace) {
            self.pos += 1;
        }
    }

    /// The next non-space character, consumed if it is one of `ops`.
    fn next_op(&mut self, ops: &[char]) -> Option<char> {
        self.skip_spaces();
        let c = self.peek().filter(|c| ops.contains(c))?;
        self.pos += 1;
        Some(c)
    }

    fn error(&self, msg: impl Into<String>) -> ExprError {
        self.error_at(self.pos, msg)
    }

    fn error_at(&self, pos: usize, msg: impl Into<String>) -> ExprError {
        ExprError {
            pos: pos + 1,
            msg: msg.into(),
        }
    }

    /// Checks the size of `value`, computed by the operator at `pos`.
    fn sized(&self, value: Integer, pos: usize) -> Result<Integer, ExprError> {
        if value.significant_bits() > MAX_BITS {
            return Err(self.too_large(pos));
        }
        Ok(value)
    }

    fn too_large(&self, pos: usize) -> ExprError {
        self.error_at(pos, format!("value too large (more than {MAX_BITS} bits)"))
    }

    /// `expr := term (('+' | '-') term)*`
    fn expr(&mut self) -> Result<Integer, ExprError> {
        let mut value = self.term()?;
        while let Some(op) = self.next_op(&['+', '-']) {
            let pos = self.pos - 1;
            let rhs = self.term()?;
            value = self.sized(if op == '+' { value + rhs } else { value - rhs }, pos)?;
        }
        Ok(value)
    }

    /// `term := unary (('*' | '/') unary)*`
    fn term(&mut self) -> Result<Integer, ExprError> {
        let mut value = self.unary()?;
        while let Some(op) = self.next_op(&['*', '/']) {
            let pos = self.pos - 1;
            let rhs = self.unary()?;
            if op == '*' {
                if value.significant_bits() + rhs.significant_bits() > MAX_BITS + 1 {
                    return Err(self.too_large(pos));
                }
                value = self.sized(value * rhs, pos)?;
            } else {
                if rhs == 0 {
                    return Err(self.error_at(pos, "division by zero"));
                }
                let (q, r) = value.div_rem_ref(&rhs).complete();
                if r != 0 {
                    return Err(self.error_at(pos, "inexact division"));
                }
                value = q;
            }
        }
        Ok(value)
    }

    /// `unary := ('-' | '+') unary | power`
    fn unary(&mut self) -> Result<Integer, ExprError> {
        if self.depth == MAX_DEPTH {
            self.skip_spaces();
            return Err(self.error(format!("nested too deeply (more than {MAX_DEPTH} levels)")));
        }
        self.depth += 1;
        let value = match self.next_op(&['-', '+']) {
            Some('-') => self.unary().map(|v| -v),
            Some(_) => self.unary(),
            None => self.power(),
        };
        self.depth -= 1;
        value
    }

    /// `power := postfix ('^' unary)?`: right associative, `-2^2 = -4`, `2^-1` is an error.
    fn power(&mut self) -> Result<Integer, ExprError> {
        let base = self.postfix()?;
        let Some(_) = self.next_op(&['^']) else {
            return Ok(base);
        };
        let pos = self.pos - 1;
        let exp_pos = self.pos;
        let exp = self.unary()?;
        if exp < 0 {
            return Err(self.error_at(exp_pos, "negative exponent"));
        }
        let bits = u64::from(base.significant_bits().saturating_sub(1));
        let fits = exp.to_u32().filter(|&e| {
            base.significant_bits() <= 1 || bits * u64::from(e) <= u64::from(MAX_BITS)
        });
        if base == 1 || base == 0 || base == -1 {
            let odd = exp.is_odd();
            return Ok(if base == -1 && !odd || base == 0 && exp == 0 {
                Integer::from(1)
            } else {
                base
            });
        }
        let Some(exp) = fits else {
            return Err(self.too_large(pos));
        };
        self.sized(base.pow(exp), pos)
    }

    /// `postfix := primary ('!' digits? | '#' digits?)*`
    fn postfix(&mut self) -> Result<Integer, ExprError> {
        let mut value = self.primary()?;
        while let Some(op) = self.next_op(&['!', '#']) {
            let pos = self.pos - 1;
            // `n!k` and `n#k`: digits right after the operator.
            let k = if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                let k_pos = self.pos;
                let k = self.number()?;
                Some(k.to_u32().filter(|&k| k > 0).ok_or_else(|| {
                    self.error_at(
                        k_pos,
                        format!("invalid step of '{op}' (must be 1 to 2^32-1)"),
                    )
                })?)
            } else {
                None
            };
            let max = if op == '!' {
                MAX_FACTORIAL
            } else {
                MAX_PRIMORIAL
            };
            let Some(n) = value.to_u32().filter(|&n| n <= max) else {
                return Err(self.error_at(
                    pos,
                    format!("argument of '{op}' must be between 0 and {max}"),
                ));
            };
            value = match (op, k) {
                ('!', None) => Integer::factorial(n).complete(),
                ('!', Some(k)) => Integer::factorial_m(n, k).complete(),
                (_, None) => Integer::primorial(n).complete(),
                (_, Some(k)) => {
                    // Product of the primes in [k, n].
                    let below = Integer::primorial(k.saturating_sub(1).min(n)).complete();
                    Integer::primorial(n).complete() / below
                }
            };
        }
        Ok(value)
    }

    /// `primary := integer | '(' expr ')'`
    fn primary(&mut self) -> Result<Integer, ExprError> {
        self.skip_spaces();
        match self.peek() {
            Some('(') => {
                let open = self.pos;
                self.pos += 1;
                let value = self.expr()?;
                if self.next_op(&[')']).is_none() {
                    self.skip_spaces();
                    return Err(match self.peek() {
                        None => self.error_at(open, "unbalanced '('"),
                        Some(c) => self.error(format!("expected ')', found '{c}'")),
                    });
                }
                Ok(value)
            }
            Some(c) if c.is_ascii_digit() => self.number(),
            Some(c) => Err(self.error(format!("expected a number, found '{c}'"))),
            None => Err(self.error("expected a number, found the end")),
        }
    }

    /// A decimal integer.
    fn number(&mut self) -> Result<Integer, ExprError> {
        let start = self.pos;
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            self.pos += 1;
        }
        if self.pos - start > (MAX_BITS / 3) as usize {
            return Err(self.too_large(start));
        }
        let digits: String = self.chars[start..self.pos].iter().collect();
        Ok(digits.parse().expect("decimal digits"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ok(input: &str) -> String {
        eval(input).unwrap().to_string()
    }

    fn err(input: &str) -> String {
        eval(input).unwrap_err().to_string()
    }

    #[test]
    fn values() {
        assert_eq!(ok("123456789"), "123456789");
        assert_eq!(ok(" 2^67 - 1 "), "147573952589676412927");
        assert_eq!(ok("2^3^2"), "512");
        assert_eq!(ok("-2^2+10"), "6");
        assert_eq!(ok("3*5+(2+7)^10"), "3486784416");
        assert_eq!(ok("(10^20-1)/9"), "11111111111111111111");
        assert_eq!(ok("10!+1"), "3628801");
        assert_eq!(ok("15!3"), "29160");
        assert_eq!(ok("11#"), "2310");
        assert_eq!(ok("17#5"), "85085");
        assert_eq!(ok("30#+1"), "6469693231");
        assert_eq!(ok("3!!"), "720");
        assert_eq!(ok("0!"), "1");
        assert_eq!(ok("1^100000000000"), "1");
        assert_eq!(ok("(0-1)^3"), "-1");
        assert_eq!(ok("2*-3"), "-6");
    }

    #[test]
    fn errors() {
        assert_eq!(err(""), "empty expression (at character 1)");
        assert_eq!(err("12a"), "unexpected 'a' (at character 3)");
        assert_eq!(err("(1+2"), "unbalanced '(' (at character 1)");
        assert_eq!(err("1+2)"), "unbalanced ')' (at character 4)");
        assert_eq!(err("7/2"), "inexact division (at character 2)");
        assert_eq!(err("7/0"), "division by zero (at character 2)");
        assert_eq!(err("2^-1"), "negative exponent (at character 3)");
        assert_eq!(
            err("1+"),
            "expected a number, found the end (at character 3)"
        );
        assert!(err("10^10^10").starts_with("value too large"));
        assert!(err("2^1048577").starts_with("value too large"));
        assert!(err("100000!").starts_with("argument of '!'"));
        assert!(err("(2^1000000)*(2^1000000)").starts_with("value too large"));
        assert!(err("1000000#").starts_with("argument of '#'"));
        assert!(err("5!0").starts_with("invalid step"));
        // Recursion is bounded: no stack overflow.
        let deep = format!("{}1{}", "(".repeat(100_000), ")".repeat(100_000));
        assert!(err(&deep).starts_with("nested too deeply"));
        assert!(err(&format!("{}5", "-".repeat(100_000))).starts_with("nested too deeply"));
        assert!(err(&vec!["2"; 100_000].join("^")).starts_with("nested too deeply"));
        let ok_depth = format!("{}7{}", "(".repeat(150), ")".repeat(150));
        assert_eq!(ok(&format!("-{ok_depth}^1")), "-7");
    }
}
