"""Codes obtained by removing several fermions from the 16.

For k = 1..4 removed fermions with independent charge vectors, counts the removed sets by the
Z-distance of the resulting [[16-k, k, d]] code (qubits use the sign bits with the all-ones row).
Checks: k = 2 always gives d = 2; for odd p, k = 3 gives d = 2 exactly for the triples inside a Weyl
image of the Pati-Salam quartet {nu^c, u^c_1, u^c_2, u^c_3}; k = 3 gives d = 1 for qubits; k = 4 gives
d = 1 for qutrits but some [[12,4,2]]_p codes for p >= 5; and nu^c - u^c_1 - u^c_2 - u^c_3 = 2 e^c.
"""
import itertools

import numpy as np

from sm16 import spinor, rank_mod, z_distance, FIELDS, NU_C, weyl_group

EXPECTED = {
    3: {1: {3: 16}, 2: {2: 120}, 3: {1: 400, 2: 160}, 4: {1: 1720}},
    5: {1: {3: 16}, 2: {2: 120}, 3: {1: 400, 2: 160}, 4: {1: 1640, 2: 80}},
    2: {1: {3: 16}, 2: {2: 120}, 3: {1: 560}, 4: {1: 1680}},
}


def table(cols, p, ks=range(1, 5)):
    out = {}
    for k in ks:
        cnt = {}
        for S in itertools.combinations(range(len(cols)), k):
            if rank_mod([cols[s] for s in S], p) < k:
                continue
            d = z_distance(cols, S, p, maxd=3) or 3
            cnt[d] = cnt.get(d, 0) + 1
        out[k] = dict(sorted(cnt.items()))
    return out


def main():
    W = spinor(5)
    f = {n: np.array(q) for n, q in FIELDS.items()}
    assert (np.array(NU_C) - f['u^c_1'] - f['u^c_2'] - f['u^c_3'] == 2 * f['e^c']).all()
    bits = [tuple([1] + [(1 - x) // 2 for x in q]) for q in W]
    for p in (3, 5, 7, 11, 2):
        t = table(bits if p == 2 else W, p)
        print(f'p={p}:', t)
        assert t == EXPECTED[p if p in EXPECTED else 5]
    quartet = [NU_C] + [FIELDS[f'u^c_{i}'] for i in (1, 2, 3)]
    images = {frozenset(g(x) for x in quartet) for g in weyl_group()}
    in_quartet = {frozenset(T) for Q4 in images for T in itertools.combinations(Q4, 3)}
    for p in (3, 5):
        good = {frozenset(W[s] for s in S) for S in itertools.combinations(range(16), 3)
                if z_distance(W, S, p, maxd=3) == 2}
        assert good == in_quartet and len(good) == 160
    print('nu^c - u^c_1 - u^c_2 - u^c_3 = 2 e^c; the 160 triples with d = 2 (p = 3, 5) are exactly '
          'the triples inside Weyl images of the quartet {nu^c, u^c_1, u^c_2, u^c_3}')


if __name__ == '__main__':
    main()
