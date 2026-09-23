#!/usr/bin/env bash
# Instructions per curve from the gungraun JSON summaries (`--save-summary=json`) of the
# `ops::curve::one_curve` benchmarks (one complete curve at the bounds of the `success_rate`
# example, benchmark id `p<digits>`), with the base branch value when there was a baseline:
# `[{"digits": 15, "pr": 8090053, "base": 8123456 | null}, ...]`.
#
# Usage: gungraun-per-curve.sh [target/gungraun] > per_curve.json
set -euo pipefail

mapfile -d '' files < <(find "${1:-target/gungraun}" -path '*/one_curve.*' -name summary.json -print0)
if [ "${#files[@]}" -eq 0 ]; then
  echo '[]'
  exit
fi
jq -s '
  def num: .Int // .Float;
  map(.profiles[0].summaries.total.summary.Callgrind.Ir.metrics as $m
    | select($m.Both or $m.Left)
    | {
        digits: (.id | ltrimstr("p") | tonumber),
        pr: (($m.Both[0] // $m.Left) | num),
        base: (if $m.Both then ($m.Both[1] | num) else null end)
      })
  | sort_by(.digits)' "${files[@]}"
