#!/usr/bin/env bash
# Per-curve wall-clock time of ecm-rs and GMP-ECM, side by side: stage 1 and stage 2 of the same
# curve (`-param 2`, same sigma) on the same numbers (balanced semiprimes of 256 and 1024 bits),
# with B1/B2 from 50k/12.7M to 1M/1e9 (see examples/per_curve.rs).
#
# Usage: scripts/compare_gmp_ecm.sh [path/to/ecm] [per_curve args, e.g. --reps 1 --bits 256]
#
# The GMP-ECM binary is the first argument, or $ECM, or `ecm` in the PATH; without one, only
# prints a message. Env: PROFILE (cargo profile of ecm-rs, default `release`; `bench-lto` for
# the optimizations of a final build). Pin both to one core for stabler numbers, e.g.
# `taskset -c 2 scripts/compare_gmp_ecm.sh`. GMP-ECM times have a 1 ms resolution, and it runs
# once per row.
set -euo pipefail

ecm=${ECM:-ecm}
if [ $# -gt 0 ] && [[ $1 != --* ]]; then
  ecm=$1
  shift
fi
if ! command -v "$ecm" > /dev/null; then
  echo "GMP-ECM not found ('$ecm'): pass its path as the first argument or in \$ECM." >&2
  exit 0
fi
# Absolute path: the script runs from the repository root.
ecm=$(realpath "$(command -v "$ecm")")

cd "$(dirname "$0")/.."
profile=${PROFILE:-release}
cargo build --quiet --profile "$profile" --features bench --example per_curve
dir=$profile
[ "$profile" = dev ] && dir=debug
per_curve=${CARGO_TARGET_DIR:-target}/$dir/examples/per_curve

ours=$(mktemp)
"$per_curve" "$@" | grep -v '^#' > "$ours"

echo "| bits | B1 / B2 | ecm-rs stage 1 + 2 (ms) | plan | GMP-ECM stage 1 + 2 (ms) | ecm-rs / GMP-ECM |"
echo "| ---: | --- | ---: | --- | ---: | ---: |"
"$per_curve" --numbers "$@" | while read -r bits b1 b2 sigma n; do
  read -r _ _ _ plan s1 s2 < <(awk -v b="$bits" -v b1="$b1" '$1 == b && $2 == b1' "$ours")
  out=$(echo "$n" | "$ecm" -c 1 -param 2 -sigma "2:$sigma" "$b1" "$b2" 2>&1 || true)
  g1=$(sed -n 's/^Step 1 took \([0-9]*\)ms.*/\1/p' <<< "$out")
  g2=$(sed -n 's/^Step 2 took \([0-9]*\)ms.*/\1/p' <<< "$out")
  ratio=$(awk -v a="$s1" -v b="$s2" -v c="${g1:-0}" -v d="${g2:-0}" \
    'BEGIN { if (c + d > 0) printf "%.2f", (a + b) / (c + d); else print "-" }')
  echo "| $bits | $b1 / $b2 | $s1 + $s2 | $plan | ${g1:-?} + ${g2:-?} | $ratio |"
done
rm -f "$ours"
