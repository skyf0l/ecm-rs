#!/usr/bin/env bash
# Converts the gungraun JSON summaries (`--save-summary=json`) into the `customSmallerIsBetter`
# format of github-action-benchmark: one entry per benchmark, with its instruction count.
#
# Usage: gungraun-to-benchmark.sh [target/gungraun] > instructions.json
set -euo pipefail

find "${1:-target/gungraun}" -name summary.json -print0 |
  xargs -0 jq -s '
    map({
      name: "\(.module_path)::\(.id)",
      unit: "instructions",
      # Callgrind instruction count (schema 7, gungraun 0.20): {"new": 123, "old": 456}.
      value: first(.profiles[] | select(.tool == "Callgrind")).data.total.metrics.Ir.values.new
    })
    | sort_by(.name)'
