// 01_generation_code.m -- the main theorem on the Standard Model code (all parts except
// uniqueness and automorphisms), the dictionary proposition, the exact-gate proposition,
// the discussion of which fermion is removed, and the spinor family.
load "common.m";
print "=== 01: one generation as a [[15,1,3]]_p code ===";

// --- The 16 weights and the anomaly ---------------------------------------
Report(Set(All16) eq Set(SpinorWeights(5)) and #All16 eq 16,
      "the 15 field columns + nu^c are exactly the 16 weights of the spinor of so(10)");
Report(#[q : q in Cols | Multiplet(q) eq "10"] eq 10 and #[q : q in Cols | Multiplet(q) eq "5bar"] eq 5,
      "15 fields split as 10 + 5bar under su(5)");
cubic := &and[ &+[Z| q[i]*q[j]*q[k] : q in All16] eq 0 : i, j, k in [1..5] ];
Report(cubic, "all cubic moments sum_q q_i q_j q_k vanish (no cubic anomaly): lattice is triorthogonal");
Report(&and[ &+[Z| q[i] : q in All16] eq 0 : i in [1..5] ], "linear moments vanish");
Report(&and[ &+[Z| q[i]*q[j] : q in All16] eq (i eq j select 16 else 0) : i, j in [1..5] ],
      "quadratic moments are 16*delta_ij");
Report(&+Nu eq 5 and
      { &+[Z| x : x in q] : q in Cols | Multiplet(q) eq "10" } eq {1} and
      { &+[Z| x : x in q] : q in Cols | Multiplet(q) eq "5bar" } eq {-3},
      "X = q.1 equals 5, 1, -3 on 1, 10, 5bar");

// Logical eigenvalue: t = e1, restricted to the 15 fields.
lam := &+[Z| q[1]^3 : q in Cols];
Report(lam eq -1, "lambda = sum_{l != nu} (e1.q_l)^3 = -1  (logical U^{-1})");
Report(&+[Z| (&+q)^3 : q in Cols] eq -125, "U(1)_chi anomaly on the 15 fields is -125 = -chi(nu^c)^3");

// Stabilizers = Cartan of su(5): t with t.1 = 0.  Over Z, tau vanishes on the full lattice,
// so check tau on the 15 columns for stabilizer basis rows and logical row.
StabT := [[1,-1,0,0,0],[0,1,-1,0,0],[0,0,1,-1,0],[0,0,0,1,-1]];
row := func< t | [ &+[Z| t[i]*q[i] : i in [1..5]] : q in Cols ] >;
Rows := [row([1,0,0,0,0])] cat [row(t) : t in StabT];       // logical first
c15 := [1 : l in [1..15]];
okT := true;
for a in [1..5] do for b in [a..5] do for cc in [b..5] do
  v := TauZ(Rows[a], Rows[b], Rows[cc], c15);
  if [a,b,cc] eq [1,1,1] then okT and:= (v eq -1); else okT and:= (v eq 0); end if;
end for; end for; end for;
Report(okT, "over Z: tau = 0 on all row triples except (L,L,L) = -1 (holds mod every p)");

// --- Distances for odd primes ------------------------------------------------
primes := [p : p in [3..101] | IsPrime(p)];
okd := true;
for p in primes do
  G := Matrix(GF(p), [[q[i] : q in Cols] : i in [1..5]]);
  okd and:= Rank(G) eq 5 and ZDistModP(Cols, Nu, p, 3) eq 3;
end for;
Report(okd, "d_Z = 3 via eq.(dist) for every odd prime p <= 101");

// Proof for all odd p: no 1 or 2 columns span nu over Q, and minors of +-1 matrices.
D2 := { Determinant(Matrix(Z, 2, 2, [a,b,c,d])) : a, b, c, d in [-1,1] };
D3 := {};
for i in [0..511] do
  bits := Intseq(i, 2); bits cat:= [0 : k in [1..9-#bits]];
  Include(~D3, Determinant(Matrix(Z, 3, 3, [1-2*x : x in bits])));
end for;
Report(D2 join D3 subset {0,2,-2,4,-4}, "2x2 and 3x3 minors of +-1 matrices lie in {0,+-2,+-4}");
nuQ := Vector(Q, Nu);
okQ := &and[ Rank(Matrix(Q, [Cols[a], Nu])) eq 2 : a in [1..15] ] and
       &and[ Rank(Matrix(Q, [Cols[a], Cols[b], Nu])) eq 3 : a, b in [1..15] | a lt b ];
Report(okQ, "over Q no one or two fields span nu  => d >= 3 for EVERY odd prime");

// --- Exact-gate proposition: exact qutrit gate, phases mod 9 ---------------
okE := true;
for t in VectorSpace(GF(3), 5) do
  ti := ChangeUniverse(Eltseq(t), Z);
  x := [ (&+[Z| ti[i]*q[i] : i in [1..5]]) mod 3 : q in All16 ];   // codeword, entries 0..2
  a := x[16];
  okE and:= ((&+[Z| x[l]^3 : l in [1..15]]) mod 9) eq ((-a^3) mod 9);
end for;
Report(okE, "every basis state of |a_bar> gets phase zeta^{-a^3}: U^{(x)15} = logical U^{-1} exactly (p=3)");


// --- Which fermion is removed ---------------------------------------------------
SMgens := func< q | [ q[1]-q[2], q[2]-q[3], (q[4]-q[5]) div 2,
                      q[1]+q[2]+q[3] - (3*(q[4]+q[5])) div 2 ] >;   // last one = 6Y
AllNames := Names cat ["nuc"];
neutralZ := [AllNames[l] : l in [1..16] | SMgens(All16[l]) eq [0,0,0,0]];
neutral3 := [AllNames[l] : l in [1..16] | [x mod 3 : x in SMgens(All16[l])] eq [0,0,0,0]];
printf "    SM-neutral over Z: %o;  mod 3: %o\n", neutralZ, neutral3;
okp := &and[ [AllNames[l] : l in [1..16] | [x mod p : x in SMgens(All16[l])] eq [0,0,0,0]] eq ["nuc"]
             : p in [p : p in [5..101] | IsPrime(p)] ];
Report(neutralZ eq ["nuc"] and Set(neutral3) eq {"ec","nuc"} and okp,
      "SM Cartan generators are stabilizers only for nu^c (Z, p>=5); {nu^c, e^c} for p=3");

// Transitivity: removing any weight gives distance 3 (p=3, spot check of Weyl transitivity).
Report(&and[ ZDistModP([q : q in All16 | q ne w], w, 3, 3) eq 3 : w in All16 ],
      "removing any one of the 16 weights gives d_Z = 3 (p=3)");

// --- Spinor of so(14): [[63,1,3]]_p -----------------------------------------------
W7 := SpinorWeights(7); nu7 := [1 : i in [1..7]];
Report(&and[ &+[Z| q[i]*q[j]*q[k] : q in W7] eq 0 : i, j, k in [1..7] ], "spinor of so(14) is anomaly free");
C63 := [q : q in W7 | q ne nu7];
Report(&and[ ZDistModP(C63, nu7, p, 3) eq 3 : p in [3,5] ], "so(14) spinor minus nu: d_Z = 3 for p = 3,5");
print "";
