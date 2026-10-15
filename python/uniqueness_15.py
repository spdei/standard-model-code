"""Uniqueness of the [[15,1,3]]_3 triorthogonal code.

Runs the exhaustive search of qutrit_triorth_search.py for every (r, m) with
r + m = 15 allowed by r + 1 <= m <= (3^r - 1)/2 - r, and checks: there are 36
double-lex solutions, all with r = 4; each has d_Z = 3; each is equivalent, by a
permutation of the qutrits and a change of basis mapping stabilizers to
stabilizers, to the code obtained by removing nu^c from the 16; and no solution
admits a second logical row (the vectors v with tau(a,b,v) = tau(a,h,v) =
tau(h,h,v) = 0 for all stabilizers a, b span exactly the stabilizer space).
Takes about ten minutes. Requires qutrit_triorth_search.py on the path.
"""
import itertools
import sys
import time

import numpy as np
import sympy

import qutrit_triorth_search as qs
from sm16 import spinor, NU_C

# The Standard Model code: logical t = e_1, stabilizers e_i - e_1
Q15 = np.array([q for q in spinor(5) if q != NU_C]).T
H_SM = np.vstack([Q15[0], Q15[1:] - Q15[0]]) % 3


def full(q, P):
    r, _ = P.shape
    return np.vstack([np.r_[np.zeros(r, int), q], np.hstack([np.eye(r, dtype=int), P])]) % 3


def invariants(H):
    """Number of weight-3 Z-logical supports through each qutrit."""
    n, A = H.shape[1], H[1:]
    c = np.zeros(n, int)
    for J in itertools.combinations(range(n), 3):
        for z in itertools.product([1, 2], repeat=3):
            v = np.zeros(n, int)
            v[list(J)] = z
            if not ((A @ v) % 3).any() and (H[0] @ v) % 3:
                c[list(J)] += 1
                break
    return c


def equivalent(Ha, Hb):
    """Find G (invertible, stabilizers to stabilizers) and a permutation with G Ha = Hb[:, pi]."""
    ca, cb = invariants(Ha), invariants(Hb)
    if sorted(ca) != sorted(cb):
        return False
    n, k = Ha.shape[1], Ha.shape[0]
    J = next(c for c in itertools.combinations(range(n), k)
             if sympy.Matrix(Ha[:, c]).det() % 3)
    Binv = np.array(sympy.Matrix(Ha[:, J]).inv_mod(3)).astype(int)
    target = sorted(tuple(Hb[:, j]) for j in range(n))
    cand = [[j for j in range(n) if cb[j] == ca[i]] for i in J]
    for imgs in itertools.product(*cand):
        if len(set(imgs)) < k:
            continue
        G = (Hb[:, imgs] @ Binv) % 3
        if (G[1:, 0] % 3).any():
            continue
        M = (G @ Ha) % 3
        if sorted(tuple(M[:, j]) for j in range(n)) == target:
            return True
    return False


def no_second_logical(H):
    """Solutions v of the linear conditions on a second logical row equal the stabilizer space."""
    h, A = H[0], H[1:]
    rows = [a * b for a, b in itertools.combinations_with_replacement(A, 2)]
    rows += [a * h for a in A] + [h * h]
    R = np.array(rows) % 3
    null_dim = H.shape[1] - qs.rank_mod3(R)
    # stabilizers satisfy the conditions; equality iff the solution space has dimension dim A
    return null_dim == qs.rank_mod3(A) and not ((R @ A.T) % 3).any()


# Search-tree node counts at n = 15 (as listed in magma/README.md).
EXPECTED_NODES = {(4, 11): 872969, (5, 10): 213694, (6, 9): 27673, (7, 8): 4272}


def main():
    sols, ranks, nodes = [], [], {}
    for r in range(1, 8):
        m = 15 - r
        if not (r + 1 <= m <= (3 ** r - 1) // 2 - r):
            continue
        t = time.time()
        S = qs.Searcher(r, m)
        found = S.run()
        print(f"(r,m)=({r},{m}): nodes {S.nodes}, solutions {len(found)}, {time.time() - t:.1f}s",
              flush=True)
        sols += [full(q, P) for q, P in found]
        ranks += [r] * len(found)
        nodes[(r, m)] = S.nodes
    ok = all(equivalent(H, H_SM) for H in sols)
    print(f"{len(sols)} double-lex solutions at n = 15; all equivalent to the SM code: {ok}")
    ok2 = all(no_second_logical(H) for H in sols)
    print(f"no normal form admits a second logical row: {ok2}")
    dists = {qs.distance_Z(H[0], H[1:], wmax=4) for H in sols}
    print(f"distances of the solutions: {dists}")
    assert ok and ok2 and len(sols) == 36 and set(ranks) == {4} and dists == {3}
    assert nodes == EXPECTED_NODES


if __name__ == "__main__":
    sys.exit(main())
