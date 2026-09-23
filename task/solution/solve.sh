#!/bin/bash
# Oracle entry point. The derivation lives in solve.py alongside this script;
# this wrapper only locates and runs it.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOLVER="$SCRIPT_DIR/solve.py"
[ -f "$SOLVER" ] || SOLVER=/solution/solve.py

python3 "$SOLVER"
