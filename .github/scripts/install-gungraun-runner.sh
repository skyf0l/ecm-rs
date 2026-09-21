#!/usr/bin/env bash
# Installs the gungraun-runner version matching the gungraun version locked in Cargo.lock of the
# current checkout, if the installed one differs: runner and library versions must be equal, and
# the base branch of a PR may lock another version.
set -euo pipefail

version=$(cargo metadata --locked --format-version 1 \
  | jq -r '.packages[] | select(.name == "gungraun") | .version')
installed=$(gungraun-runner --version 2>/dev/null | awk '{ print $2 }' || true)

if [ "$installed" != "$version" ]; then
  echo "Installing gungraun-runner $version (installed: ${installed:-none})"
  cargo install gungraun-runner --version "$version" --locked --force
fi
