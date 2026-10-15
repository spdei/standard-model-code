"""Distillation with the fifteen-qutrit and fifteen-qubit codes: the figures of the paper's
distillation table and text.

Noise model: each input suffers a nontrivial Z^k with total probability eps, all p-1 powers equally
likely. A round accepts when the Z-error lies in A^perp and fails when it lies in A^perp minus
H^perp. Depolarizing strength delta gives eps = (p-1)/p * delta.
"""
import itertools

import numpy as np

from sm16 import spinor, NU_C


def nullspace_vectors(M, p):
    """All vectors v in F_p^n with M v = 0 (M given as a list of rows)."""
    M = [list(map(lambda x: int(x) % p, r)) for r in M]
    rows, n = len(M), len(M[0])
    piv, r = [], 0
    for c in range(n):
        k = next((i for i in range(r, rows) if M[i][c]), None)
        if k is None:
            continue
        M[r], M[k] = M[k], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    basis = []
    for f in (c for c in range(n) if c not in piv):
        v = [0] * n
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-M[i][f]) % p
        basis.append(v)
    Bm = np.array(basis, dtype=np.int64)
    for coef in itertools.product(range(p), repeat=len(basis)):
        yield (np.array(coef) @ Bm) % p


def enumerators(A, H, p):
    """Hamming weight enumerators of accepted errors (A^perp) and logical errors (A^perp minus H^perp)."""
    n = len(A[0])
    acc, log = np.zeros(n + 1, dtype=np.int64), np.zeros(n + 1, dtype=np.int64)
    Hm = np.array(H) % p
    for v in nullspace_vectors(A, p):
        w = int(np.count_nonzero(v))
        acc[w] += 1
        if ((Hm @ v) % p).any():
            log[w] += 1
    return acc, log


def round_map(acc, log, eps, p):
    n = len(acc) - 1
    q = eps / (p - 1)
    P = sum(acc[w] * (1 - eps) ** (n - w) * q ** w for w in range(n + 1))
    F = sum(log[w] * (1 - eps) ** (n - w) * q ** w for w in range(n + 1))
    return F / P, P


def threshold(acc, log, p):
    lo, hi = 1e-6, 0.999 * (p - 1) / p
    for _ in range(200):
        mid = (lo + hi) / 2
        if round_map(acc, log, mid, p)[0] < mid:
            lo = mid
        else:
            hi = mid
    return lo


def codes():
    Q = np.array([q for q in spinor(5) if q != NU_C]).T
    H3 = [list(Q[0])] + [list(Q[i] - Q[0]) for i in range(1, 5)]   # logical q_1, stabilizers
    B = (1 - Q) // 2                                                # qubits: sign bits
    A2 = [list(B[i]) for i in range(5)]
    return (H3[1:], H3, 3), (A2, [[1] * 15] + A2, 2)


# Values printed in the paper: low-weight enumerators (accepted, logical), weight-3 supports,
# leading coefficients at equal eps and equal delta, thresholds, acceptance and two-round output.
PAPER = {
    3: dict(acc=[1, 0, 0, 50, 260], log3=50, supports=25, lead_eps='6.25', lead_delta='2.78',
            eps_star='0.230', delta_star='0.345', acc_eps='0.860', acc_delta='0.905', two='1.7e-15'),
    2: dict(acc=[1, 0, 0, 35, 105], log3=35, supports=35, lead_eps='35', lead_delta='8.75',
            eps_star='0.141', delta_star='0.283', acc_eps='0.860', acc_delta='0.928', two='1.6e-12'),
}

if __name__ == '__main__':
    for (A, H, p), name in zip(codes(), ['[[15,1,3]]_3', '[[15,1,3]]_2']):
        acc, log = enumerators(A, H, p)
        lead = log[3] / (p - 1) ** 3                       # eps_out ~ lead * eps^3
        lead_d = lead * ((p - 1) / p) ** 3 * p / (p - 1)   # delta_out ~ lead_d * delta^3
        e = threshold(acc, log, p)
        x = 0.01
        for _ in range(2):
            x = round_map(acc, log, x, p)[0]
        got = dict(acc=[int(v) for v in acc[:5]], log3=int(log[3]), supports=int(log[3]) // (p - 1),
                   lead_eps=f'{lead:.4g}', lead_delta=f'{lead_d:.3g}', eps_star=f'{e:.3f}',
                   delta_star=f'{e * p / (p - 1):.3f}',
                   acc_eps=f'{round_map(acc, log, 0.01, p)[1]:.3f}',
                   acc_delta=f'{round_map(acc, log, 0.01 * (p - 1) / p, p)[1]:.3f}', two=f'{x:.2g}')
        print(name, got)
        assert got == PAPER[p], (got, PAPER[p])
    print('distillation figures match the paper')
