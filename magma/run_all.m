// run_all.m -- run from this directory:  magma run_all.m
// Expected: every line [ok]. Script 05 is the slow one (a few minutes; the n=15, r=4 search
// visits 872,969 nodes, matching the table of search pairs in README.md).
load "01_generation_code.m";
load "02_qubits.m";
load "03_dirac_doubling.m";
load "04_logicals_distillation.m";
load "05_minimality_uniqueness.m";
load "06_builtin_automorphisms.m";
quit;
