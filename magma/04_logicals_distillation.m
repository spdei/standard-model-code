// 04_logicals_distillation.m -- the minimal logical errors, the paper's table of the 25
// minimum-weight logicals as four-fermion relations,
// and the distillation comparison (weight enumerators, threshold, acceptance).
load "common.m";
print "=== 04: minimal logicals and distillation ===";
K := GF(3);

// --- 25 integer relations nu^c = sum eps_f q(f), three fields, eps = +-1 -----------
rels := [];
for J in Subsets({1..15}, 3) do
  Js := Setseq(J);
  for e in [[a,b,c] : a, b, c in [-1,1]] do
    if &and[ &+[Z| e[i]*Cols[Js[i]][k] : i in [1..3]] eq 1 : k in [1..5] ] then
      Append(~rels, <Js, e>);
    end if;
  end for;
end for;
Report(#rels eq 25, "exactly 25 relations nu^c = +-f +-f' +-f'' over Z");

// Weight-3 Z-logical supports mod 3 (no pairs, since d = 3).
V := [Vector(K, q) : q in Cols]; s := Vector(K, Nu);
sup3 := { J : J in Subsets({1..15}, 3) |
          Rank(Matrix([V[j] : j in J])) eq Rank(VerticalJoin(Matrix([V[j] : j in J]), Matrix([s]))) };
Report(#sup3 eq 25 and sup3 eq { Set(R[1]) : R in rels },
      "the 25 relations are exactly the 25 weight-3 Z-logical supports (p=3)");

// Interaction structure.
tenw := func< v | Sort([Abs(x) : x in v]) eq [0,0,0,0,2] >;          // weight of the 10 (scaled)
root := func< v | Sort([Abs(x) : x in v]) eq [0,0,0,2,2] >;          // root of so(10) (scaled)
n10 := 0; ng := 0; nmix := 0; agree := true;
print "    relations  nu^c = sum eps f :";
for R in rels do
  Js := R[1]; e := R[2];
  ten := #[ i : i in [1..3] | e[i] eq -1 and tenw([Nu[k] + Cols[Js[i]][k] : k in [1..5]]) ] gt 0;
  gb  := #[ i : i in [1..3] | e[i] eq  1 and root([Cols[Js[i]][k] - Nu[k] : k in [1..5]]) ] gt 0;
  mix := #Set(e) eq 2;
  n10 +:= ten select 1 else 0; ng +:= gb select 1 else 0; nmix +:= mix select 1 else 0;
  agree and:= (gb eq mix);
  printf "      %o%o\n", &cat[ (e[i] eq 1 select " +" else " -") cat Names[Js[i]] : i in [1..3] ],
         gb select "   [g]" else "";
end for;
Report(n10 eq 25, "every relation pairs nu^c with a field into a weight of the 10 (Yukawa exchange)");
Report(ng eq 15 and nmix eq 15 and agree, "15 relations admit gauge-boson exchange = exactly the mixed-sign ones");

// --- Weight enumerators and distillation ----------------------------------------------
G := Matrix(K, [[q[i] : q in Cols] : i in [1..5]]);
A := Matrix(K, [[ &+[Z| t[i]*q[i] : i in [1..5]] : q in Cols ] :
                t in [[1,-1,0,0,0],[0,1,-1,0,0],[0,0,1,-1,0],[0,0,0,1,-1]] ]);
CA := LinearCode(A); CH := LinearCode(G);
WA := [0 : i in [0..15]]; WH := [0 : i in [0..15]];
for t in WeightDistribution(Dual(CA)) do WA[t[1]+1] := t[2]; end for;
for t in WeightDistribution(Dual(CH)) do WH[t[1]+1] := t[2]; end for;
printf "    A^perp weight distribution: %o\n    H^perp weight distribution: %o\n", WA, WH;
Report(WA[2] eq 0 and WA[3] eq 0 and WA[4] - WH[4] eq 50 and WH[4] eq 0,
      "50 weight-3 logical Z-errors (25 supports x 2 multiples): eps_out ~ 50 (eps/2)^3 = 6.25 eps^3");

RR := RealField(30);
Pr := func< W, e | &+[ RR!W[k+1]*(1-e)^(15-k)*(e/2)^k : k in [0..15] ] >;
Acc := func< e | Pr(WA, e) >;
Out := func< e | (Pr(WA, e) - Pr(WH, e))/Pr(WA, e) >;
e0 := RR!1/100;
printf "    at eps = 0.01: acceptance = %o, eps_out = %o, eps_out/eps^3 = %o\n",
       RealField(6)!Acc(e0), RealField(6)!Out(e0), RealField(6)!(Out(e0)/e0^3);
lo := RR!1/10; hi := RR!45/100;
for i in [1..80] do
  mid := (lo + hi)/2;
  if Out(mid) lt mid then lo := mid; else hi := mid; end if;
end for;
printf "    threshold eps* = %o\n", RealField(6)!lo;
Report(Abs(lo - 0.230) lt 0.0005 and Abs(Acc(e0) - 0.86) lt 0.005, "threshold 0.230 and acceptance 0.86 at eps = 0.01");
print "";
