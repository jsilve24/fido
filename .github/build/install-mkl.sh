#!/usr/bin/env bash
set -euo pipefail

if [[ $(uname -s) != Linux || $(uname -m) != x86_64 || $# != 1 ]]; then
  echo 'Usage (Linux x86_64 only): bash install-mkl.sh ABSOLUTE_PREFIX' >&2
  exit 1
fi
prefix=$1
if [[ $prefix != /* || $prefix == *[[:space:]]* ]]; then
  echo 'The CI MKL prefix must be absolute and contain no whitespace.' >&2
  exit 1
fi
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python3 -m venv "$prefix"
"$prefix/bin/python" -m pip install --disable-pip-version-check --no-deps \
  --only-binary=:all: --require-hashes -r "$script_dir/mkl-requirements.txt"
test -f "$prefix/include/mkl.h"
test -f "$prefix/lib/libmkl_intel_lp64.so.2"
test -f "$prefix/lib/libmkl_gnu_thread.so.2"
test -f "$prefix/lib/libmkl_core.so.2"
echo "Installed CI MKL 2024.2.2 (LP64, GNU threading) at $prefix"
