"""Shared definitions: spinor weights, the field labels of the 16, the Weyl group of so(10),
and linear algebra modulo p."""
import itertools

import numpy as np


def spinor(m):
    """Weights of a spinor of so(2m), scaled by 2: sign vectors with an even number of minus signs."""
    return [q for q in itertools.product([1, -1], repeat=m) if q.count(-1) % 2 == 0]


NU_C = (1, 1, 1, 1, 1)

# The fifteen fields in the order of the paper's table of the sixteen weights (columns l and l+7
# are Dirac partners, the last is the seed nu). Coordinates 1-3 are colour, 4-5 weak.
FIELDS = {
    'd^c_1': (1, -1, -1, -1, -1), 'd^c_2': (-1, 1, -1, -1, -1), 'u^c_3': (-1, -1, 1, 1, 1),
    'd^c_3': (-1, -1, 1, -1, -1), 'u^c_2': (-1, 1, -1, 1, 1), 'u^c_1': (1, -1, -1, 1, 1),
    'e^c': (1, 1, 1, -1, -1),
    'd_1': (-1, 1, 1, -1, 1), 'd_2': (1, -1, 1, -1, 1), 'u_3': (1, 1, -1, 1, -1),
    'd_3': (1, 1, -1, -1, 1), 'u_2': (1, -1, 1, 1, -1), 'u_1': (-1, 1, 1, 1, -1),
    'e': (-1, -1, -1, -1, 1), 'nu': (-1, -1, -1, 1, -1),
}
NAMES = list(FIELDS)


def weyl_group():
    """The Weyl group of so(10): coordinate permutations and even numbers of sign changes,
    as functions on weights (1920 elements)."""
    signs = spinor(5)
    return [lambda q, P=P, s=s: tuple(s[i] * q[P[i]] for i in range(5))
            for P in itertools.permutations(range(5)) for s in signs]


def rank_mod(M, p):
    M = [list(map(lambda x: int(x) % p, r)) for r in M]
    if not M:
        return 0
    r, rows, cols = 0, len(M), len(M[0])
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
    return r


def z_distance(cols, removed, p, maxd=4):
    """Smallest |J| such that span(h_J) meets span(h_removed) nontrivially (the paper's distance
    formula, extended to several removed columns)."""
    S = [cols[s] for s in removed]
    k = len(S)
    rest = [l for l in range(len(cols)) if l not in removed]
    for d in range(1, maxd + 1):
        for J in itertools.combinations(rest, d):
            hJ = [cols[j] for j in J]
            rJ = rank_mod(hJ, p)
            if rJ == d and rank_mod(hJ + S, p) < rJ + k:
                return d
    return None


def tau_vanishes(H, c=None, logical=0, lam=None):
    """True if the integer trilinear form with coefficients c vanishes on all row triples of H,
    except (logical, logical, logical), which must equal lam (if lam is given)."""
    H = np.array(H)
    c = np.ones(H.shape[1], int) if c is None else np.array(c)
    for a, b, d in itertools.combinations_with_replacement(range(len(H)), 3):
        t = int((c * H[a] * H[b] * H[d]).sum())
        if (a, b, d) == (logical,) * 3:
            if lam is not None and t != lam:
                return False
        elif t:
            return False
    return True


def primes(lo, hi):
    return [p for p in range(lo, hi + 1) if p > 1 and all(p % k for k in range(2, int(p ** .5) + 1))]
