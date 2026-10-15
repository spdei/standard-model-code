#!/bin/sh
# Runs every Python check (about 20 minutes, mostly uniqueness_15.py). Any failure stops the run.
set -e
cd "$(dirname "$0")"
python3 letter_checks.py
python3 doubling_checks.py
python3 interactions_15.py
python3 removals_16.py
python3 spinor_checks.py
python3 distillation.py
python3 automorphisms_15.py
python3 automorphisms_p.py
python3 qutrit_triorth_search.py
python3 uniqueness_15.py
echo "all checks passed"
