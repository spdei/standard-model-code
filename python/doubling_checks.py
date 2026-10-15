"""Doubling as Dirac pairing: the doubling theorem and its proof, including the printed
seven-qudit P-code and doubled generator matrix, and the change of basis to the printed
generator matrix of the main theorem.
"""
import itertools
from fractions import Fraction as Fr

import numpy as np

from sm16 import NU_C, FIELDS, NAMES, tau_vanishes

F = [np.array(FIELDS[n]) for n in NAMES]                     # columns 0-6 f^c, 7-13 f, 14 nu
nuc = np.array(NU_C)
sigma = lambda q: np.array([-q[0], -q[1], -q[2], q[3], -q[4]])

# --- sigma and the Pati-Salam sectors ------------------------------------------------------------
assert [-1, -1, -1, 1, -1].count(-1) % 2 == 0               # sigma is in the Weyl group
assert all((sigma(F[l]) == F[l + 7]).all() for l in range(7)) and (sigma(nuc) == F[14]).all()
Qem = lambda q: Fr(q[0] + q[1] + q[2], 6) - Fr(q[4], 2)      # T_3L + Y
BL = lambda q: Fr(q[0] + q[1] + q[2], 3)
assert all(Qem(F[l]) == -Qem(F[l + 7]) and BL(F[l]) == -BL(F[l + 7]) for l in range(7))
R = F[:7] + [nuc]
assert all(q[3] == q[4] for q in R) and all(q[3] == -q[4] for q in F[7:])   # T_3L = 0 / T_3R = 0
print('sigma (a Weyl element) maps column l to l+7 and nu^c to nu; Dirac partners have opposite '
      'Q_em and B-L; T_3L vanishes on the right-handed sector and T_3R on its image')

# --- The right-handed P-code -------------------------------------------------------------------------
EQ_STEANE = np.array([[1, 0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 0, 1, 1], [0, 0, 0, 1, 1, 1, 1]])
MU = np.array([-1, -1, 1, -1, 1, 1, -1])
Rm = np.array(F[:7]).T
assert (MU == Rm[3]).all()                                     # mu = q_4 = 2 T_3R
S = np.array([(1 - Rm[3] * Rm[i]) // 2 for i in range(3)])
assert (S == EQ_STEANE).all()
assert sorted(int(S[0, l] + 2 * S[1, l] + 4 * S[2, l]) for l in range(7)) == list(range(1, 8))
assert not (S @ MU).any() and not np.einsum('il,jl,l->ij', S, S, MU).any()
assert MU.sum() == -1
for p in (3, 5, 7):                                            # P-code distance 3
    d = min(k for k in (1, 2, 3) for J in itertools.combinations(range(7), k)
            for v in itertools.product(range(1, p), repeat=k)
            if not (S[:, list(J)] @ np.array(v) % p).any() and sum(v) % p)
    assert d == 3
v = np.array([1, 1, -1, 0, 0, 0, 0])
assert not (S @ v).any() and v.sum() == 1
assert ((v * MU) @ np.array(F[:7]) == nuc).all()               # undo rescaling: nu^c relation
print('P-code: stabilizers and mu = 2T_3R as printed, columns are 1..7 in binary, P-code '
      'conditions over Z, N = -1, distance 3 (p = 3, 5, 7); logical (1,1,-1,0,...) is '
      'nu^c = -d^c_1 - d^c_2 - u^c_3')

# --- The doubled generator matrix ---------------------------------------------------------------------
one, z = np.ones(7, int), np.zeros(7, int)
EQ_DOUBLED = np.vstack([np.r_[-one, one, 1]] + [np.r_[s, -s, 0] for s in EQ_STEANE]
                       + [np.r_[z, one, 1]])
q = np.array(F).T
rows = np.vstack([-q[4]] + [(q[4] - q[i]) // 2 for i in range(3)] + [(q[3] - q[4]) // 2])
r = np.r_[q[3, :7], q[3, 7:14], 1]                             # rescale by q_4, seed by 1
assert (rows * r == EQ_DOUBLED).all()
assert tau_vanishes(EQ_DOUBLED, c=np.r_[MU, MU, 1], lam=1)
assert MU.sum() * 1 ** 3 + 1 == 0                               # N kappa^3 + lambda_T = 0
print('rescaled charges -q_5, (q_5 - q_i)/2, 2T_3L give the printed doubled matrix over Z; it is '
      'triorthogonal with coefficients mu, mu, 1 and logical exponent 1; N + lambda_T = 0')

# --- Change of basis to the main theorem's matrix --------------------------------------------------
e = np.eye(5, dtype=object)
half = Fr(1, 2)
Lt, Kt = -e[4], (e[3] - e[4]) * half
st = [(e[4] - e[i]) * half for i in range(3)]
target = [e[0]] + [(e[i] - e[i + 1]) * half for i in range(4)]
M = [[-1, -2, 0, 0, 0], [0, -1, 1, 0, 0], [0, 0, -1, 1, 0], [0, 0, 0, -1, -1], [0, 0, 0, 0, 1]]
basis = [Lt] + st + [Kt]
for row, t in zip(M, target):
    assert all(x == y for x, y in zip(sum(c * b for c, b in zip(row, basis)), t))
assert round(abs(np.linalg.det(np.array(M, float)))) == 1
assert all(x == 0 for x in (e[0] + Lt) + 2 * st[0])           # q_1 = -L modulo stabilizers
assert all(x == y for x, y in zip(-sum(st) / 3, [Fr(1, 6)] * 3 + [0, -half]))   # Q_em
print('change of basis to the printed main-theorem matrix is unimodular; q_1 = -L modulo '
      'stabilizers; Q_em = -(1/3) sum_i (q_5 - q_i)/2')

# --- Qubits: the doubled Steane code ----------------------------------------------------------------
B = (1 - np.array(F).T) // 2
gens = np.vstack([np.ones(15, int), B])
span = {tuple((np.array(c) @ gens) % 2) for c in itertools.product(range(2), repeat=6)}
assert tuple([1] * 15) in span and tuple([0] * 7 + [1] * 8) in span
ss = [w for w in span if w[:7] == w[7:14] and w[14] == 0 and any(w)]
assert len(ss) == 7 and {sum(w[:7]) for w in ss} == {4}
print('p = 2: in the same columns the qubit code contains 1, the linking vector (0,1,1) and (s,s,0) '
      'for s in the [7,3,4] simplex code (doubled Steane code)')
