// 03_dirac_doubling.m -- the doubling theorem (doubling = Dirac pairing), the paper's table
// of the sixteen weights as fields (with the doubled-form matrix H15), the P-code
// conditions, and integral doubling.
load "common.m";
print "=== 03: doubling is Dirac pairing ===";

// sigma : q -> (-q1,-q2,-q3,q4,-q5) maps f^c (column l) to its partner f (column l+7).
sigma := func< q | [-q[1], -q[2], -q[3], q[4], -q[5]] >;
Report(&and[ sigma(Cols[l]) eq Cols[l+7] : l in [1..7] ], "sigma maps column l to column l+7 (Dirac pairs)");
Report(sigma(Nu) eq Cols[15], "sigma maps nu^c to the seed nu");

// Pati-Salam charges (as rationals) and the SM charges on the 16 left-handed fields.
T3L := func< q | (q[4]-q[5])/4 >; T3R := func< q | (q[4]+q[5])/4 >;
BL  := func< q | (q[1]+q[2]+q[3])/3 >;
Y   := func< q | BL(q)/2 - T3R(q) >; Qem := func< q | T3L(q) + Y(q) >;
Report(&and[ T3L(Cols[l]) eq 0 : l in [1..7] ] and T3L(Nu) eq 0, "T3L vanishes on right-handed fields");
Report(Qem(Nu) eq 0 and Y(Nu) eq 0 and Y(Cols[7]) eq 1 and Qem(Cols[14]) eq -1,
      "Y(nu^c) = 0, Y(e^c) = 1, Q(e) = -1");
Report(&and[ Qem(Cols[l]) eq -Qem(Cols[l+7]) and BL(Cols[l]) eq -BL(Cols[l+7]) : l in [1..7] ],
      "Dirac partners have opposite Q_em and B-L");

// --- Rescaled matrix: r = q4 on copy 1, -q5 on copy 2, 1 on the seed ---------------
r := [ Cols[l][4] : l in [1..7] ] cat [ -Cols[l][5] : l in [8..14] ] cat [1];
Lrow := [ -Cols[l][5]*r[l] : l in [1..15] ];
Srows := [ [ (Cols[l][i] - Cols[l][5])*r[l] : l in [1..15] ] : i in [1..3] ];
Krow := [ ((Cols[l][4] - Cols[l][5]) div 2)*r[l] : l in [1..15] ];
Report(Lrow eq [-1 : l in [1..7]] cat [1 : l in [1..8]], "logical row -q5 becomes (-1, 1, 1)");
Report(Krow eq [0 : l in [1..7]] cat [1 : l in [1..8]], "linking row 2T3L becomes (0, 1, 1)");
Report(&and[ S[l+7] eq -S[l] : l in [1..7], S in Srows ] and &and[ S[15] eq 0 : S in Srows ],
      "rows q_i - q5 become (s, -s, 0)");
mu := r[1..7];
Report(mu eq [ Z!(2*T3R(Cols[l])) : l in [1..7] ] and r[8..14] eq mu,
      "weights after rescaling: mu = 2*T3R on both copies (+1 on u^c, -1 on d^c, e^c)");
s := [ S[1..7] : S in Srows ];
okP := &and[ &+[Z| mu[l]*s[i][l] : l in [1..7]] eq 0 : i in [1..3] ] and
       &and[ &+[Z| mu[l]*s[i][l]*s[j][l] : l in [1..7]] eq 0 : i, j in [1..3] ];
Report(okP and &+mu eq -1, "P-code conditions hold over Z, N = sum mu = -1 = -2T3R(nu^c)");

// Doubled lattice is triorthogonal over Z with weights r, lambda = 1.
Rws := [Lrow] cat Srows cat [Krow];
okD := true;
for a in [1..5] do for b in [a..5] do for c in [b..5] do
  v := TauZ(Rws[a], Rws[b], Rws[c], r);
  okD and:= (([a,b,c] eq [1,1,1]) select (v eq 1) else (v eq 0));
end for; end for; end for;
Report(okD, "doubled matrix (kappa=1) is an integral triorthogonal lattice, tau(L,L,L) = 1");

// P-code distance d_S = min{wt v : v in S^perp, sum v != 0} over GF(3).
K := GF(3); SM := Matrix(K, s);
dS := Min([ Weight(v) : v in VectorSpace(K, 7) | IsZero(SM*Transpose(Matrix([v]))) and &+Eltseq(v) ne 0 ]);
Report(dS eq 3, "P-code [[7,1,3]]^P_3 has distance 3; doubling distance min(3, 1+2) = 3");

// --- H15 rows and comparison with the Standard Model code ---------------
H15c := [ [0,0,0,0,-1], [0,0,1,0,-1], [0,1,1,0,-2], [1,0,1,0,-2], [0,0,0,-1,1] ];
ev := func< c, q | &+[Z| c[i]*q[i] : i in [1..5]] >;   // -q5, q3-q5, q2+q3-2q5, q1+q3-2q5, q5-q4
H15 := Matrix(K, [ [ ev(c, q) : q in Cols ] : c in H15c ]);
print "    H15 (mod 3), columns in the H15 column order:"; print H15;
Gthm := Matrix(K, [ [ q[i] : q in Cols ] : i in [1..5] ]);
Athm := Matrix(K, [ [ ev(t, q) : q in Cols ] :
                    t in [[1,-1,0,0,0],[0,1,-1,0,0],[0,0,1,-1,0],[0,0,0,1,-1]] ]);
Report(RowSpace(H15) eq RowSpace(Gthm), "row space of H15 = row space of the Standard Model code");
Report(&and[ ev(c, Nu) eq 0 : c in H15c[2..5] ] and ev(H15c[1], Nu) eq -1 and
      RowSpace(RowSubmatrix(H15, 2, 4)) eq RowSpace(Athm),
      "H15 rows 2-5 span the su(5) Cartan stabilizers; row 1 is logical (value -1 on nu^c)");

// --- Over Q: Cartan of SU(3)xU(1)_em, logical row = -(B-L), B-L = chi/5 mod stabilizers ---
// Functionals on charge space as coefficient vectors c, acting by q -> c.q.
VQ := VectorSpace(Q, 5);
span := sub< VQ | [VQ![1,0,0,0,-1], VQ![0,1,0,0,-1], VQ![0,0,1,0,-1]] >;
Report(VQ![1/6,1/6,1/6,0,-1/2] in span and VQ![1,-1,0,0,0] in span and VQ![0,1,-1,0,0] in span,
      "q_i - q5 span Q_em and the colour Cartan (over Q)");
nuv := VQ!Nu;
isStab := func< c | InnerProduct(c, nuv) eq 0 >;
Report(isStab(VQ![0,0,0,0,-1] - (-VQ![1/3,1/3,1/3,0,0])), "-q5 = -(B-L) modulo stabilizers");
Report(isStab(VQ![1/3,1/3,1/3,0,0] - VQ![1/5,1/5,1/5,1/5,1/5]), "B-L = chi/5 modulo stabilizers");
Report(isStab(VQ![1,0,0,0,0] - VQ![1/5,1/5,1/5,1/5,1/5]), "t = e1 is chi/5 modulo stabilizers");
print "";
