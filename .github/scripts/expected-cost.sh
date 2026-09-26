#!/usr/bin/env bash
# Markdown table of the expected cost to find a factor of each size, for the PR comment and the
# job summary: expected curves (`success_rate` example, `--json`) x instructions per curve
# (`one_curve` benchmarks, see gungraun-per-curve.sh), for the PR and, when both of its inputs
# exist, the base branch. A change is only marked better or worse when the 95% intervals
# (those of the expected curves, times the instructions per curve) don't overlap.
#
# Usage: expected-cost.sh per_curve.json success_rate.json [base_success_rate.json] > cost.md
set -euo pipefail

# Missing inputs (e.g. no baseline on the base branch) count as no results.
empty=$(mktemp)
echo '[]' > "$empty"
input() { if [ -n "${1:-}" ] && [ -f "$1" ]; then echo "$1"; else echo "$empty"; fi; }
per_curve=$(input "${1:-}")
success=$(input "${2:-}")
base_success=$(input "${3:-}")

jq -n -r --slurpfile ir "$per_curve" --slurpfile pr "$success" --slurpfile base "$base_success" '
  def human:
    if . >= 1e12 then "\(. / 1e12 * 100 | round / 100)T"
    elif . >= 1e9 then "\(. / 1e9 * 100 | round / 100)G"
    elif . >= 1e6 then "\(. / 1e6 * 100 | round / 100)M"
    elif . >= 1e3 then "\(. / 1e3 * 100 | round / 100)K"
    else tostring end;
  def one: . * 10 | round / 10;
  # {digits, bounds, curves, low, high} of a success rate result.
  def curves:
    (.name | capture("(?<d>[0-9]+)-digit factor \\((?<b>[^)]*)\\)")) as $n
    | (.extra | capture("95% CI: (?<lo>[0-9.]+)\\.\\.(?<hi>[0-9.]+|inf)")) as $ci
    | {digits: ($n.d | tonumber), bounds: $n.b, curves: .value, low: ($ci.lo | tonumber),
       high: (if $ci.hi == "inf" then infinite else ($ci.hi | tonumber) end)};
  # Expected cost of `c` (curves) with `ir` instructions per curve, or null.
  def cost(c; ir):
    if c == null or ir == null then null
    else {value: (c.curves * ir), low: (c.low * ir), high: (c.high * ir), curves: c.curves, ir: ir}
    end;
  def cell: if . == null then "-"
    else "\(.curves | one) x \(.ir | human) = **\(.value | human)**"
      + " (\(.low | human)..\(if .high == infinite then "inf" else (.high | human) end))" end;

  ($ir[0] // []) as $ir
  | [($pr[0] // [])[] | curves] as $pr
  | [($base[0] // [])[] | curves] as $base
  | [ $pr[]
      | . as $p
      | ($ir | map(select(.digits == $p.digits)) | first) as $i
      | select($i != null)
      | cost($p; $i.pr) as $new
      | cost($base | map(select(.digits == $p.digits)) | first; $i.base) as $old
      | {name: "\($p.digits) digits (\($p.bounds))", new: $new, old: $old} ] as $rows
  | "### Expected cost to find a factor\n",
    if ($rows | length) == 0 then
      "_No results: needs the `one_curve` instruction counts and the curve success rate._"
    else
      "Instructions to find a factor: expected curves (success rate) x instructions per curve"
        + " (`one_curve`, one curve at the same bounds, on a number of the same size). 95%"
        + " intervals from the expected curves.\n",
      "| | Factor | Base | PR | Change |",
      "|---|---|---:|---:|---:|",
      ($rows[]
        | (if .old == null then ""
           elif .new.high < .old.low then "🟢"
           elif .new.low > .old.high then "🔴"
           else "⚪" end) as $mark
        | (if .old == null then "new"
           else ((.new.value - .old.value) / .old.value * 100 | one) as $v
             | if $v > 0 then "+\($v)%" else "\($v)%" end end) as $change
        | "| \($mark) | \(.name) | \(.old | cell) | \(.new | cell) | \($change) |")
    end'
