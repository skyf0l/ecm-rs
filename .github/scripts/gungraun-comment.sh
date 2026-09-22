#!/usr/bin/env bash
# Markdown comparison of the gungraun results with their baseline (`--baseline`), for the PR
# comment: sorted from the worst to the best instruction count change, unchanged benchmarks
# collapsed. The regression limit, the big improvement mark and the noise level are in percent.
#
# Usage: gungraun-comment.sh [target/gungraun] [limit] [noise] [big] > instructions.md
set -euo pipefail

find "${1:-target/gungraun}" -name summary.json -print0 |
  xargs -0 jq -s -r --argjson limit "${2:-2}" --argjson noise "${3:-0.5}" --argjson big "${4:-10}" '
    def num: .Int // .Float;
    def human:
      if . >= 1e9 then "\(. / 1e9 * 100 | round / 100)G"
      elif . >= 1e6 then "\(. / 1e6 * 100 | round / 100)M"
      elif . >= 1e3 then "\(. / 1e3 * 100 | round / 100)K"
      else tostring end;
    def pct: (. * 100 | round / 100) as $v
      | if $v == 0 then "0%" elif $v > 0 then "+\($v)%" else "\($v)%" end;
    def change(m): (m.metrics.Both[0] | num) as $new | (m.metrics.Both[1] | num) as $old
      | {new: $new, old: $old, pct: (($new - $old) / $old * 100)};
    def row: "| \(.name) | \(.ir.old | human) → \(.ir.new | human) | \(.ir.pct | pct) | \(.cycles.pct | pct) |";
    def mark:
      if .ir.pct > $limit then "🔴"
      elif .ir.pct > $noise then "🟠"
      elif .ir.pct <= -$big then "🚀"
      else "🟢" end;

    # Benchmarks without a base value (added by the PR).
    ([ .[]
      | .profiles[0].summaries.total.summary.Callgrind as $c
      | select($c.Ir.metrics.Left)
      | {name: "`\(.function_name)` \(.id)", ir: ($c.Ir.metrics.Left | num)} ]
      | sort_by(.name)) as $new

    | [ .[]
      | .profiles[0].summaries.total.summary.Callgrind as $c
      | select($c.Ir.metrics.Both)
      | {name: "`\(.function_name)` \(.id)", ir: change($c.Ir), cycles: change($c.EstimatedCycles)} ]
    | sort_by(-.ir.pct)
    | map(select(.ir.pct > $noise)) as $worse
    | map(select(.ir.pct < -$noise)) as $better
    | map(select(.ir.pct | fabs <= $noise)) as $same
    | "### Instruction counts\n",
      if length == 0 then
        "_No baseline on the base branch: nothing to compare with._"
      else
        if ($worse + $better | length) > 0 then
          "| | Benchmark | Instructions (base → PR) | Change | Est. cycles |",
          "|---|---|---:|---:|---:|",
          (($worse + $better)[] | "| \(mark) " + row),
          ""
        else empty end,
        if ($new | length) > 0 then
          "🆕 \($new | length) new: " + ($new | map("\(.name) (\(.ir | human))") | join(", ")) + "\n"
        else empty end,
        if ($same | length) > 0 then
          "<details><summary>⚪ \($same | length) unchanged (within ±\($noise)%)</summary>\n",
          "| Benchmark | Instructions (base → PR) | Change | Est. cycles |",
          "|---|---:|---:|---:|",
          ($same[] | row),
          "\n</details>"
        else empty end
      end'
