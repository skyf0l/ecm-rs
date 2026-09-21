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
      value: (.profiles[0].summaries.total.summary.Callgrind.Ir.metrics | (.Both[0] // .Left).Int)
    })
    | sort_by(.name)'
