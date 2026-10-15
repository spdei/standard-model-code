"""Monomial automorphism group of the fifteen-qutrit code (p = 3): order 240, consisting of the
permutations of the five coordinates of q together with an overall sign (S_5 x Z_2), with orbits
the 10 and the 5bar.

An automorphism is a signed permutation of the qutrits mapping the stabilizer space A and the row
space H to themselves. It preserves the column invariants (number of weight-three logical supports
through each qutrit), so it is determined by the signed images of a fixed information set among
columns with equal invariants; each candidate is tested on all fifteen columns.
"""
import itertools
import numpy as np
import sympy
from sm16 import spinor, NU_C, rank_mod

W = [q for q in spinor(5) if q != NU_C]
Q = np.array(W).T
H = np.vstack([Q[0], Q[1:] - Q[0]]) % 3          # logical row q_1, stabilizers q_i - q_1


def invariants(H):
    n, A = H.shape[1], H[1:]
    c = np.zeros(n, int)
    for J in itertools.combinations(range(n), 3):
        for z in itertools.product([1, 2], repeat=3):
            v = np.zeros(n, int); v[list(J)] = z
            if not ((A @ v) % 3).any() and (H[0] @ v) % 3:
                c[list(J)] += 1
                break
    return c


def main():
    n = 15
    c = invariants(H)
    J = next(cc for cc in itertools.combinations(range(n), 5) if sympy.Matrix(H[:, cc]).det() % 3)
    Binv = np.array(sympy.Matrix(H[:, J]).inv_mod(3)).astype(int)
    target = {tuple(H[:, j]): j for j in range(n)}
    cand = [[j for j in range(n) if c[j] == c[i]] for i in J]
    count = 0
    for imgs in itertools.product(*cand):
        if len(set(imgs)) < 5:
            continue
        for sg in itertools.product([1, 2], repeat=5):
            G = ((H[:, imgs] * np.array(sg)) @ Binv) % 3
            if (G[1:, 0] % 3).any():                  # stabilizers must map to stabilizers
                continue
            M = (G @ H) % 3
            used, ok = set(), True
            for j in range(n):
                k = target.get(tuple(M[:, j]), target.get(tuple((-M[:, j]) % 3)))
                if k is None or k in used:
                    ok = False; break
                used.add(k)
            if ok:
                count += 1
    print('order of the monomial automorphism group (p = 3):', count)
    assert count == 240
    # The 120 coordinate permutations, with an overall sign, give 240 distinct automorphisms.
    maps = set()
    for P in itertools.permutations(range(5)):
        perm = [W.index(tuple(q[P[i]] for i in range(5))) for q in W]
        for s in (1, 2):
            M = (H[:, perm] * s) % 3
            assert rank_mod(np.vstack([H, M]).tolist(), 3) == 5
            assert rank_mod(np.vstack([H[1:], M[1:]]).tolist(), 3) == 4
            maps.add((tuple(perm), s))
    assert len(maps) == 240
    orbits = {frozenset(W.index(tuple(q[P[i]] for i in range(5))) for P in itertools.permutations(range(5)))
              for q in W}
    assert sorted(len(o) for o in orbits) == [5, 10]
    print('these are the coordinate permutations with an overall sign: S_5 x Z_2, orbits 10 and 5bar')


if __name__ == '__main__':
    main()
