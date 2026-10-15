"""One generation and the [[15,1,3]]_p code: the claims of the main text up to and including the
main theorem, the qubit case, the choice of removed fermion, and the single-U(1) example.
"""
import itertools
from fractions import Fraction as Fr

import numpy as np

from sm16 import spinor, NU_C, FIELDS, NAMES, z_distance, tau_vanishes, primes

W16 = spinor(5)
Q15 = [FIELDS[n] for n in NAMES]
assert sorted(Q15 + [NU_C]) == sorted(W16)

# --- The 16 and its anomalies ---------------------------------------------------------------------
A16 = np.array(W16)
assert not np.einsum('li,lj,lk->ijk', A16, A16, A16).any()              # cubic sums
assert not A16.sum(axis=0).any()                                        # linear sums
assert (A16.T @ A16 == 16 * np.eye(5, dtype=int)).all()                 # quadratic sums
assert all(tuple(-x for x in q) not in W16 for q in W16)                # chiral
assert {q.count(-1): sum(q) for q in W16} == {0: 5, 2: 1, 4: -3}        # chi on 1, 10, 5bar
assert sum(sum(q) ** 3 for q in Q15) == -125 == -sum(NU_C) ** 3        # U(1)_chi on the 15
print('16: cubic and linear sums vanish, quadratic sums 16 delta_ij, chiral; chi = 5, 1, -3; '
      'chi anomaly of the 15 is -125')

# --- Field labels and charges (table of the sixteen weights) ----------------------------------------
def charges(q):
    nc = sum(q[i] == -1 for i in range(3))
    nw = sum(q[i] == -1 for i in (3, 4))
    Y = Fr(nw, 2) - Fr(nc, 3)
    T3L, T3R, BL = Fr(q[3] - q[4], 4), Fr(q[3] + q[4], 4), Fr(q[0] + q[1] + q[2], 3)
    return Y, T3L, T3R, BL

SM_Y = {'u^c': Fr(-2, 3), 'd^c': Fr(1, 3), 'e^c': 1, 'u': Fr(1, 6), 'd': Fr(1, 6), 'e': Fr(-1, 2),
        'nu': Fr(-1, 2)}
SM_T3L = {'u': Fr(1, 2), 'd': Fr(-1, 2), 'nu': Fr(1, 2), 'e': Fr(-1, 2)}
for name, q in FIELDS.items():
    base = name.split('_')[0] if name[0] in 'ud' else name
    Y, T3L, T3R, BL = charges(q)
    assert Y == SM_Y[base] and T3L == SM_T3L.get(base, 0), name
    assert Y == BL / 2 - T3R                                            # Y = (B-L)/2 - T_3R
    if name[0] in 'ud':                                                 # colour subscript
        mark = 1 if '^c' in name else -1                                # antiquark: unoccupied
        assert [i for i in range(3) if q[i] == mark] == [int(name[-1]) - 1], name
assert charges(NU_C)[0] == 0
print('field labels: Y, T_3L, colour and Y = (B-L)/2 - T_3R agree with the Standard Model')

# --- The generator matrix printed in the paper ------------------------------------------------------
EQ_H = np.array([
    [-1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1],
    [0, -1, -1, -1, 1, 1, 1, 0, 0, 0, 1, -1, 0, 0, 0],
    [-1, 1, 0, 0, 0, -1, -1, 1, 1, 0, 0, 1, -1, 0, 0],
    [0, -1, 1, 0, -1, 1, 0, 0, -1, 1, 0, 0, 1, -1, 0],
    [0, 0, -1, 1, 0, -1, 1, -1, 1, 0, 0, 0, 0, 1, -1]])
EQ_H_FIELDS = ['u^c_3', 'u^c_2', 'd_1', 'u_1', 'u^c_1', 'd_2', 'u_2', 'd_3', 'u_3', 'e^c',
               'd^c_1', 'd^c_2', 'd^c_3', 'nu', 'e']
cols = [tuple(-1 if k in (i, j) else 1 for k in range(5))              # 10: minus signs at ij
        for i, j in itertools.combinations(range(5), 2)]
cols += [tuple(1 if t == k else -1 for t in range(5)) for k in range(5)]  # 5bar: plus sign at k
assert cols == [FIELDS[n] for n in EQ_H_FIELDS]
C = np.array(cols).T
assert (EQ_H == np.vstack([C[0]] + [(C[i] - C[i + 1]) // 2 for i in range(4)])).all()
H16 = np.hstack([np.array([[1], [0], [0], [0], [0]]), EQ_H])            # restore nu^c
assert tau_vanishes(H16)
coef = np.array([np.linalg.lstsq(EQ_H.T.astype(float), C[i].astype(float), rcond=None)[0]
                 for i in range(5)])
assert np.allclose(coef, np.round(coef)) and np.allclose(np.round(coef) @ EQ_H, C)
assert round(abs(np.linalg.det(np.round(coef)))) == 16                  # index 2^4
print('printed generator matrix: equals q_1, (q_i - q_{i+1})/2 on the stated columns and labels; '
      'triorthogonal over Z with nu^c restored; contains the charge lattice with index 16')

# --- The code for odd p ------------------------------------------------------------------------------
A = EQ_H[1:]
assert tau_vanishes(EQ_H, lam=-1)                                       # lambda = -1 over Z
for p in primes(3, 101):
    assert z_distance(W16, (W16.index(NU_C),), p) == 3, p
for k in (1, 2):                                                        # over Q
    for J in itertools.combinations(Q15, k):
        assert np.linalg.matrix_rank(np.array(list(J) + [NU_C])) == k + 1
for n, expected in ((2, {0, 2, -2}), (3, {0, 4, -4})):
    assert {round(np.linalg.det(np.array(M).reshape(n, n)))
            for M in itertools.product([1, -1], repeat=n * n)} == expected
print('odd p: triorthogonal with lambda = -1; d_Z = 3 for all odd p <= 101; over Q no one or two '
      'weights span nu^c; 2x2 and 3x3 minors of +-1 matrices lie in {0, +-2} and {0, +-4}')

def exact(H, p, mod):
    """Sum over the fifteen qudits of x^3 equals -a^3 modulo mod, for every codeword."""
    rep = lambda x: x if x <= p // 2 else x - p
    for c in itertools.product(range(p), repeat=5):
        x = [rep(int(v)) for v in (np.array(c) @ H) % p]
        if (sum(v ** 3 for v in x[1:]) + x[0] ** 3) % mod:
            return False
    return True
assert exact(H16, 3, 9) and exact(H16, 5, 5) and exact(H16, 7, 7)
assert int(A[0] @ A[0]) == 8
print('transversal gate is logical U^dagger exactly (modulo 9 for p = 3; p = 5, 7); '
      'not self-orthogonal for any odd p (a stabilizer row has norm 8)')

# --- Logical operator and the removed fermion -------------------------------------------------------
BL_t = [Fr(1, 3)] * 3 + [0, 0]
chi5_t = [Fr(1, 5)] * 5
assert sum(BL_t) == sum(chi5_t) == 1                                    # charge 1 at nu^c
SM_gens = lambda q: [q[0] - q[1], q[1] - q[2], (q[3] - q[4]) // 2,
                     q[0] + q[1] + q[2] - 3 * (q[3] + q[4]) // 2]       # colour, 2T_3L, 6Y
assert all(SM_gens(q)[3] == 6 * charges(q)[0] for q in W16)
ALL = list(FIELDS.items()) + [('nu^c', NU_C)]
assert [n for n, q in ALL if not any(SM_gens(q))] == ['nu^c']
assert sorted(n for n, q in ALL if all(g % 3 == 0 for g in SM_gens(q))) == ['e^c', 'nu^c']
for p in primes(5, 101):
    assert [n for n, q in ALL if all(g % p == 0 for g in SM_gens(q))] == ['nu^c']
for w in W16:                                                           # Weyl transitivity
    assert sorted(tuple(a * b for a, b in zip(w, q)) for q in W16) == sorted(W16)
print('B-L and chi/5 differ by a stabilizer over Q; the Standard Model Cartan generators are all '
      'stabilizers only for nu^c (over Z and mod 5 <= p <= 101), and also for e^c mod 3; '
      'the Weyl group is transitive on the 16')

f = {**{n: np.array(q) for n, q in FIELDS.items()}, 'nu^c': np.array(NU_C)}
assert (f['nu^c'] == f['e^c'] + f['u^c_1'] - f['d^c_1']).all()
assert (f['nu^c'] == f['u_2'] + f['u^c_2'] - f['nu']).all()
assert (f['nu^c'] == -f['d^c_1'] - f['d^c_2'] - f['u^c_3']).all()
print('relations nu^c = e^c + u^c_1 - d^c_1, u_2 + u^c_2 - nu and -d^c_1 - d^c_2 - u^c_3 hold')

# --- Qubits ------------------------------------------------------------------------------------------
B = (1 - A16.T) // 2                                                    # sign bits, 5 x 16
n = 16
S1, S2, S3 = A16.sum(0), A16.T @ A16, np.einsum('li,lj,lk->ijk', A16, A16, A16)
for i, j, k in itertools.permutations(range(5), 3):
    assert B[i].sum() == 8 == (n - S1[i]) // 2
    assert (B[i] * B[j]).sum() == 4 == (n - S1[i] - S1[j] + S2[i, j]) // 4
    assert (B[i] * B[j] * B[k]).sum() == 2 == (n - S1[i] - S1[j] - S1[k] + S2[i, j] + S2[i, k]
                                               + S2[j, k] - S3[i, j, k]) // 8
gens = np.vstack([np.ones(16, int), B])
span = {tuple((np.array(c) @ gens) % 2) for c in itertools.product(range(2), repeat=6)}
x4 = [tuple((1 - np.array(q[:4])) // 2) for q in W16]                  # points of F_2^4
rm14 = {tuple((a + np.dot(c, x)) % 2 for x in x4)
        for a in range(2) for c in itertools.product(range(2), repeat=4)}
assert len(set(x4)) == 16 and span == rm14 and {sum(v) for v in span} == {0, 8, 16}
s0 = W16.index(NU_C)
assert not B[:, s0].any()
keep = [l for l in range(16) if l != s0]
assert all(sum(v[l] for l in keep) % 8 == (-v[s0]) % 8 for v in span)  # T^15 = logical T^dagger
A2 = B[:, keep]
dz2 = min(len(J) for k in (1, 2, 3) for J in itertools.combinations(range(15), k)
          if not (A2[:, list(J)].sum(1) % 2).any() and len(J) % 2)
assert dz2 == 3
print('qubits: sign-bit weights 8, overlaps 4 and 2 (directly and from the charge sums); with 1 '
      'they span RM(1,4); removing nu^c gives the 15-qubit code, d_Z = 3, T^15 = logical T^dagger')

# --- The single-U(1) example ----------------------------------------------------------------------------
x = [1, 5, -7, -8, 9]
assert sum(x) == 0 and sum(v ** 3 for v in x) == 0 and sum(v ** 3 for v in x[1:]) == -1
print('example (1, 5, -7, -8, 9): linear and cubic sums vanish; logical exponent -1')
