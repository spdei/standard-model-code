"""The 25 minimum-weight logicals as four-fermion relations nu^c = sum_f eps_f q(f), eps_f = +-1.

Checks the eight families and their counts in the paper's table of the 25 relations; that the
15 mixed-sign relations are exactly those splitting into two pairs joined by a root (gauge-boson
exchange); that every relation splits into two pairs joined by a weight of the 10, with nu^c in one
pair; which component of the 10 each can exchange; and the baryon- and lepton-number statements.
"""
import itertools
from collections import Counter
from fractions import Fraction as Fr

import numpy as np

from sm16 import FIELDS, NU_C

NUC = np.array(NU_C)
Q = {k: np.array(v) for k, v in FIELDS.items()}
E = np.eye(5, dtype=int)
TEN = [s * 2 * E[i] for i in range(5) for s in (1, -1)]
ROOTS = [2 * (a * E[i] + b * E[j]) for i, j in itertools.combinations(range(5), 2)
         for a in (1, -1) for b in (1, -1)]
isin = lambda v, L: any((v == w).all() for w in L)
base = lambda f: f.split('_')[0] if f[0] in 'ud' else f
colour = lambda f: int(f[-1]) if f[0] in 'ud' else None
B = lambda f: Fr(1, 3) * (1 if '^c' not in f else -1) if f[0] in 'ud' else Fr(0)
L = lambda f: {'e': 1, 'nu': 1, 'e^c': -1}.get(f, 0)

# (name, signed fields, colour pattern, count, mixed signs) -- the paper's table
FAMILIES = [
    ('(nu nu^c)(u u^c)^dagger', [(1, 'u'), (1, 'u^c'), (-1, 'nu')], 'same', 3, True),
    ('(nu nu^c)(d d^c)', [(-1, 'd'), (-1, 'd^c'), (-1, 'nu')], 'same', 3, False),
    ('(nu nu^c)(e e^c)', [(-1, 'e'), (-1, 'e^c'), (-1, 'nu')], None, 1, False),
    ('right-handed current', [(1, 'e^c'), (1, 'u^c'), (-1, 'd^c')], 'same', 3, True),
    ('mixed current', [(1, 'd'), (1, 'u^c'), (-1, 'e')], 'same', 3, True),
    ('Delta B = Delta L, d u d^c', [(1, 'd'), (1, 'u'), (-1, 'd^c')], 'distinct', 6, True),
    ('Delta B = Delta L, d^c d^c u^c', [(-1, 'd^c'), (-1, 'd^c'), (-1, 'u^c')], 'distinct', 3, False),
    ('scalar charged current', [(-1, 'u'), (-1, 'd^c'), (-1, 'e')], 'same', 3, False),
]


def relations():
    out = []
    for J in itertools.combinations(FIELDS, 3):
        for eps in itertools.product([1, -1], repeat=3):
            if (sum(e * Q[f] for e, f in zip(eps, J)) == NUC).all():
                out.append((J, eps))
    return out


def family(J, eps):
    sig = sorted((e, base(f)) for e, f in zip(eps, J))
    cols = [colour(f) for f in J if colour(f) is not None]
    pattern = None if not cols else 'same' if len(set(cols)) == 1 else \
        'distinct' if len(set(cols)) == len(cols) else 'other'
    hits = [name for name, s, pat, _, _ in FAMILIES if sorted(s) == sig and pat == pattern]
    assert len(hits) == 1, (J, eps)
    return hits[0]


def main():
    rels = relations()
    assert len(rels) == 25
    counts = Counter(family(J, eps) for J, eps in rels)
    assert counts == {name: n for name, _, _, n, _ in FAMILIES}
    for J, eps in rels:
        print(' nu^c =', ' '.join(f"{'+' if e > 0 else '-'}{f}" for e, f in zip(eps, J)),
              ' [' + family(J, eps) + ']')
        ops = [(1, NUC)] + [(-e, Q[f]) for e, f in zip(eps, J)]   # nu^c with three other fields
        mixed = len(set(eps)) > 1
        assert mixed == dict((n, m) for n, _, _, _, m in FAMILIES)[family(J, eps)]
        # every relation pairs nu^c with a field into a weight of the 10
        assert any(isin(NUC + s * q, TEN) for s, q in ops[1:])
        # gauge-boson exchange (two pairs joined by a root) exactly for the mixed-sign relations
        assert any(isin(NUC - q, ROOTS) for s, q in ops[1:] if s == -1) == mixed
        dB = sum(e * B(f) for e, f in zip(eps, J))
        dL = sum(e * L(f) for e, f in zip(eps, J)) + 1                    # L(nu^c) = -1
        assert dB == dL                                                   # B - L conserved
        fam = family(J, eps)
        assert (dB != 0) == fam.startswith('Delta B')
        if fam == 'scalar charged current':
            assert dB == dL == 0
    print('25 relations in the eight families of the table, with its counts; 15 mixed-sign, exactly '
          'those admitting gauge-boson exchange; all admit exchange of a 10; B - L conserved, '
          'B violated by the 9 Delta B = Delta L relations only')
    exchange_components(rels)


def exchange_components(rels):
    """A weak coordinate (i = 4, 5) of the exchanged 10 weight is a Higgs doublet, a colour
    coordinate (i = 1, 2, 3) a colour triplet."""
    for J, eps in rels:
        kinds = set()
        for f, e in zip(J, eps):
            v = NUC - e * Q[f]
            if isin(v, TEN):
                kinds.add('triplet' if np.nonzero(v)[0][0] < 3 else 'doublet')
        fam = family(J, eps)
        needs_triplet = fam == 'right-handed current' or fam.startswith('Delta B')
        assert (kinds == {'triplet'}) == needs_triplet, (J, eps, kinds)
    print('right-handed currents and Delta B relations require a colour triplet; the Dirac-mass, '
          'mixed and scalar relations can arise from Higgs-doublet exchange')


if __name__ == '__main__':
    main()
