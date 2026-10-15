// 02_qubits.m -- the qubit criterion and its proof: sign bits, RM(1,4),
// the 15-qubit code with transversal T, and the doubled Steane structure.
load "common.m";
print "=== 02: qubits from sign bits ===";
F2 := GF(2);
B := [ [ (1 - q[i]) div 2 : q in All16 ] : i in [1..5] ];      // sign bits, length 16

ok1 := &and[ &+B[i] eq 8 and (16 - &+[Z| q[i] : q in All16]) div 2 eq 8 : i in [1..5] ];
ok2 := &and[ &+[Z| B[i][l]*B[j][l] : l in [1..16]] eq 4 and
             (16 + &+[Z| q[i]*q[j] : q in All16]) div 4 eq 4 : i, j in [1..5] | i lt j ];
ok3 := &and[ &+[Z| B[i][l]*B[j][l]*B[k][l] : l in [1..16]] eq 2 and
             2 - (&+[Z| q[i]*q[j]*q[k] : q in All16]) div 8 eq 2 : i, j, k in [1..5] | i lt j and j lt k ];
Report(ok1 and ok2 and ok3, "|b_i| = 8, |b_i^b_j| = 4, |b_i^b_j^b_k| = 2, also via the charge-moment formulas");

C16 := LinearCode(Matrix(F2, [[1 : l in [1..16]]] cat B));
Report(Dimension(C16) eq 5 and IsEquivalent(C16, ReedMullerCode(1, 4)), "span(1, b_i) = RM(1,4)");
Report(WeightDistribution(C16) eq [<0,1>,<8,30>,<16,1>], "weights 0, 8, 16: triply even");

// Remove nu^c (index 16, where every b_i vanishes).
Report(&and[ B[i][16] eq 0 : i in [1..5] ], "all sign bits vanish at nu^c");
Stab := LinearCode(Matrix(F2, [ B[i][1..15] : i in [1..5] ]));
Hq   := LinearCode(Matrix(F2, [ B[i][1..15] : i in [1..5] ] cat [[1 : l in [1..15]]]));
Report(Dimension(Stab) eq 4 and Dimension(Hq) eq 5, "15-qubit code: 4 X-stabilizers + logical all-ones");
okT := &and[ (w in Stab) select (Weight(w) mod 8 eq 0) else (Weight(w) mod 8 eq 7) : w in Hq ];
Report(okT, "stabilizer coset weights = 0 mod 8, logical coset = 7 mod 8: T^(x)15 = logical T^dagger");
dZ := Min([ Weight(w) : w in Dual(Stab) | not (w in Dual(Hq)) ]);
Report(dZ eq 3, "Z-distance of the 15-qubit code is 3");

// Doubled Steane: columns 1-7 = right-handed fields, 8-14 = Dirac partners, 15 = nu (seed).
Report(Vector(F2, [1 : l in [1..15]]) in Hq, "contains all-ones");
Report(Vector(F2, [0 : l in [1..7]] cat [1 : l in [1..8]]) in Hq, "contains linking vector (0, 1, 1)");
Sset := [ s : s in VectorSpace(F2, 7) | Vector(F2, Eltseq(s) cat Eltseq(s) cat [0]) in Hq ];
Scode := LinearCode(sub< VectorSpace(F2, 7) | Sset >);
Report(#Sset eq 8 and Dimension(Scode) eq 3 and
      {Weight(s) : s in Sset | not IsZero(s)} eq {4},
      "{s : (s,s,0) in code} is a [7,3,4] simplex code => doubled Steane code");
print "";
