# Magma scripts

These scripts require [Magma](http://magma.maths.usyd.edu.au/magma/), which needs a licence. They are an optional cross-check of the Python scripts.

These scripts do not use the Python code. Run them from this folder:

```
cd magma
magma run_all.m
```

Each check prints `[ok]` or `[FAIL]` followed by a description.

| File | Topic |
|---|---|
| `common.m` | Shared definitions |
| `01_generation_code.m` | The 16, the code for odd p and its gate, the choice of removed fermion, the so(14) spinor |
| `02_qubits.m` | The qubit case |
| `03_dirac_doubling.m` | Doubling and Dirac pairing |
| `04_logicals_distillation.m` | Minimum-weight logicals and distillation figures |
| `05_minimality_uniqueness.m` | Exhaustive search, uniqueness, automorphism group (p = 3) |
| `06_builtin_automorphisms.m` | Automorphism groups with `AutomorphismGroup` (p = 3, 5, 7, 11) |
| `run_all.m` | Runs scripts 01 to 06 |

Search-tree node counts printed by script 05, for each (n, r, m) allowed by the bounds:

| n | r | m | nodes | solutions |
|---|---|---|---:|---:|
| 7 | 3 | 4 | 18 | 0 |
| 8 | 3 | 5 | 54 | 0 |
| 9 | 3 | 6 | 192 | 0 |
| 9 | 4 | 5 | 55 | 0 |
| 10 | 3 | 7 | 694 | 0 |
| 10 | 4 | 6 | 198 | 0 |
| 11 | 3 | 8 | 2,882 | 0 |
| 11 | 4 | 7 | 831 | 0 |
| 11 | 5 | 6 | 199 | 0 |
| 12 | 3 | 9 | 12,385 | 0 |
| 12 | 4 | 8 | 3,877 | 0 |
| 12 | 5 | 7 | 838 | 0 |
| 13 | 3 | 10 | 52,721 | 0 |
| 13 | 4 | 9 | 23,229 | 0 |
| 13 | 5 | 8 | 4,263 | 0 |
| 13 | 6 | 7 | 839 | 0 |
| 14 | 4 | 10 | 134,003 | 0 |
| 14 | 5 | 9 | 26,551 | 0 |
| 14 | 6 | 8 | 4,271 | 0 |
| 15 | 4 | 11 | 872,969 | 36 |
| 15 | 5 | 10 | 213,694 | 0 |
| 15 | 6 | 9 | 27,673 | 0 |
| 15 | 7 | 8 | 4,272 | 0 |
