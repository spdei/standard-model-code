// 05_minimality_uniqueness.m -- minimality and uniqueness:
//   * no triorthogonal qutrit code with d >= 3 and n <= 14,
//   * all [[15,1,3]]_3 codes equivalent to the Standard Model code,
//   * no [[15,k,3]]_3 with k >= 2, no [[15,1,4]]_3,
//   * monomial automorphism group S5 x Z2 with orbits 10 + 5.
// Normal form (the paper's reduction to a normal form): A = (I_r | P), h = (0 | q), unit weights,
// conditions (a)-(g).  Search over double-lex representatives: q = (0^z 1^o 2^t) with
// lambda = o - t = 1 mod 3, rows of P lex non-decreasing, columns of (q;P) lex non-decreasing.
// The search is iterative (explicit stack); node counts should reproduce the table of search pairs
// (r,m) in magma/README.md, which python/search_output.txt also reproduces for n <= 14.
load "common.m";
print "=== 05: minimality, uniqueness, automorphisms (p = 3) ===";
F3 := GF(3);

function SearchRM(r, m)
  sols := [];  nodes := 0;
  for z in [0..m] do for o in [0..m-z] do
    t := m - z - o;
    if (o - t) mod 3 ne 1 then continue; end if;
    q := [0 : i in [1..z]] cat [1 : i in [1..o]] cat [2 : i in [1..t]];
    // Candidate rows, in lex order, satisfying the single-row conditions (a), (e), (f, i=j).
    RI := [];
    for i in [0..3^m-1] do
      d := Intseq(i, 3); d := Reverse(d cat [0 : k in [1..m-#d]]);
      if (&+d) mod 3 ne 2 then continue; end if;                                  // (a)
      if (&+[Z| d[l] : l in [1..m] | q[l] ne 0]) mod 3 ne 0 then continue; end if; // (e)
      if (&+[Z| q[l]*d[l]^2 : l in [1..m]]) mod 3 ne 0 then continue; end if;      // (f) i=j
      Append(~RI, d);
    end for;
    N := #RI;
    if N eq 0 then continue; end if;
    Aug := [ Vector(F3, d cat [x^2 : x in d]) : d in RI ];     // (p | p*p)
    stC := [ [Z|] ];  stK := [ [Z| j : j in [1..N]] ];  stS := [ [ q[l] eq q[l+1] : l in [1..m-1] ] ];
    while #stC gt 0 do
      chosen := stC[#stC]; cand := stK[#stK]; same := stS[#stS];
      Prune(~stC); Prune(~stK); Prune(~stS);
      nodes +:= 1;
      if #chosen eq r then
        // (g): columns of P of weight >= 2 and pairwise non-parallel.
        okg := true; seen := {};
        for l in [1..m] do
          col := [ RI[i][l] : i in chosen ];
          nz := [x : x in col | x ne 0];
          if #nz lt 2 then okg := false; break; end if;
          f := nz[1]; nc := [ (x*f) mod 3 : x in col ];          // f^{-1} = f in F3
          if nc in seen then okg := false; break; end if;
          Include(~seen, nc);
        end for;
        if okg then Append(~sols, <q, Matrix(F3, [RI[i] : i in chosen])>); end if;
        continue;
      end if;
      for k in [1..#cand] do
        idx := cand[k]; d := RI[idx];
        okm := true;                                          // columns stay lex ordered
        for l in [1..m-1] do
          if same[l] and d[l] gt d[l+1] then okm := false; break; end if;
        end for;
        if not okm then continue; end if;
        newsame := [ same[l] and d[l] eq d[l+1] : l in [1..m-1] ];
        zero := [0 : l in [1..m]];
        cons := [ Vector(F3, [x^2 : x in d] cat zero),                    // (b): sum_{supp p} p' = 0
                  Vector(F3, [q[l]*d[l] : l in [1..m]] cat zero),         // (f): sum q p p' = 0
                  Vector(F3, zero cat d) ];                               // (b): sum_{supp p'} p = 0
        for i in chosen do
          Append(~cons, Vector(F3, [RI[i][l]*d[l] : l in [1..m]] cat zero)); // (c)
        end for;
        rest := cand[k..#cand];
        Prod := Matrix([ Aug[j] : j in rest ]) * Transpose(Matrix(cons));
        newcand := [Z| rest[j] : j in [1..#rest] | IsZero(Prod[j]) ];
        Append(~stC, chosen cat [idx]); Append(~stK, newcand); Append(~stS, newsame);
      end for;
    end while;
  end for; end for;
  return sols, nodes;
end function;

// Pairs (r,m) allowed by the bounds r+1 <= m <= (3^r-1)/2 - r, n = r+m <= 15.
pairs := [ <r, n-r> : n in [7..15], r in [3..7] | r+1 le n-r and n-r le (3^r-1) div 2 - r ];
printf "    %o pairs with n <= 14, %o with n = 15\n",
       #[P : P in pairs | P[1]+P[2] le 14], #[P : P in pairs | P[1]+P[2] eq 15];
printf "    %3o %3o %3o %9o %5o %8o\n", "n", "r", "m", "nodes", "sols", "time";
Sols15 := [];  short := 0;
for P in pairs do
  r := P[1]; m := P[2]; t0 := Cputime();
  sols, nodes := SearchRM(r, m);
  printf "    %3o %3o %3o %9o %5o %8o\n", r+m, r, m, nodes, #sols, Cputime(t0);
  if r+m le 14 then short +:= #sols; else Sols15 cat:= [<r, S> : S in sols]; end if;
end for;
Report(short eq 0, "no triorthogonal qutrit code with d >= 3 and n <= 14");
Report(#Sols15 eq 36 and { S[1] : S in Sols15 } eq {4}, "n = 15: 36 double-lex solutions, all with r = 4");

// --- Build each solution as 15 column vectors in F3^5 (coordinate 1 = logical row) ----
function SolCols(S)
  q := S[1]; P := S[2]; r := Nrows(P); m := Ncols(P);
  cols := [ Vector(F3, [0] cat [ (i eq j) select 1 else 0 : j in [1..r]]) : i in [1..r] ];
  cols cat:= [ Vector(F3, [F3!q[l]] cat [P[i][l] : i in [1..r]]) : l in [1..m] ];
  return cols;
end function;

function W3Counts(cols, s)
  n := #cols; cnt := [0 : i in [1..n]];
  for J in Subsets({1..n}, 3) do
    M := Matrix([cols[j] : j in J]);
    if Rank(M) eq Rank(VerticalJoin(M, Matrix([s]))) then
      for j in J do cnt[j] +:= 1; end for;
    end if;
  end for;
  return cnt;
end function;

// All M in GL(5,3) with cols[l]*M in tgt (or +-tgt if signed) for every l, and s*M = +-sT.
// Search over images of an information set, restricted by the invariant inv.
function FindMaps(cols, s, inv, tgt, sT, invT, signed)
  I := [Z|];
  for j in [1..#cols] do
    if Rank(Matrix([cols[i] : i in I cat [j]])) eq #I + 1 then Append(~I, j); end if;
    if #I eq 5 then break; end if;
  end for;
  Binv := Matrix([cols[i] : i in I])^-1;
  tset := {@ v : v in tgt @};
  sg := signed select [F3!1, F3!-1] else [F3!1];
  C := [ [ <e, k> : k in [1..#tgt], e in sg | invT[k] eq inv[I[i]] ] : i in [1..5] ];
  found := [];
  for c1 in C[1] do for c2 in C[2] do
    if c2[2] eq c1[2] then continue; end if;
    for c3 in C[3] do
      if c3[2] in {c1[2], c2[2]} then continue; end if;
      for c4 in C[4] do
        if c4[2] in {c1[2], c2[2], c3[2]} then continue; end if;
        for c5 in C[5] do
          if c5[2] in {c1[2], c2[2], c3[2], c4[2]} then continue; end if;
          Tm := Matrix([ c[1]*tgt[c[2]] : c in [c1, c2, c3, c4, c5] ]);
          if Rank(Tm) lt 5 then continue; end if;
          M := Binv * Tm; ok := true;
          for v in cols do
            w := v*M;
            if not (w in tset or (signed and (-w) in tset)) then ok := false; break; end if;
          end for;
          if ok and (s*M eq sT or s*M eq -sT) then Append(~found, M); end if;
        end for;
      end for;
    end for;
  end for; end for;
  return found;
end function;

Tcols := [ Vector(F3, q) : q in Cols ]; sT := Vector(F3, Nu);
invT := W3Counts(Tcols, sT);
printf "    invariants of the Standard Model code (weight-3 logicals per qutrit): %o\n", invT;

okEq := true; okD := true; okK2 := true; okTri := true;
for S in Sols15 do
  cols := SolCols(S[2]); s := Vector(F3, [1,0,0,0,0]);
  // sanity: triorthogonality with lambda = 1, and distance exactly 3
  rowsS := [ Vector(F3, [c[i] : c in cols]) : i in [1..5] ];   // row 1 = h, rows 2..5 = stabilizers
  for a in [1..5] do for b in [a..5] do for c in [b..5] do
    v := &+Eltseq(CP(CP(rowsS[a], rowsS[b]), rowsS[c]));
    okTri and:= (([a,b,c] eq [1,1,1]) select (v eq 1) else (v eq 0));
  end for; end for; end for;
  okD and:= ZDistVec(cols, s, 3) eq 3;
  // equivalence to the Standard Model code (exact, unsigned column maps)
  inv := W3Counts(cols, s);
  if Sort(inv) ne Sort(invT) then okEq := false; continue; end if;
  maps := FindMaps(cols, s, inv, Tcols, sT, invT, false);
  okEq and:= #maps gt 0;
  // k >= 2: second logical row v needs tau(a,b,v) = tau(a,h,v) = tau(h,h,v) = 0
  h := rowsS[1]; As := rowsS[2..5];
  cons := [ CP(As[i], As[j]) : i, j in [1..4] | i le j ] cat [ CP(a, h) : a in As ] cat [ CP(h, h) ];
  NS := Nullspace(Transpose(Matrix(cons)));
  okK2 and:= Dimension(NS) eq 4 and NS eq sub< VectorSpace(F3, 15) | As >;
end for;
Report(okTri, "every solution is triorthogonal with lambda = 1");
Report(okD, "every solution has d = 3 exactly (so no [[15,1,4]]_3)");
Report(okEq, "all 36 solutions are equivalent to the Standard Model code  => uniqueness");
Report(okK2, "for all 36: admissible second logical rows span only the stabilizers => no [[15,k,3]]_3, k >= 2");

// --- Monomial automorphism group of the Standard Model code ----------------------------------
auts := FindMaps(Tcols, sT, invT, Tcols, sT, invT, true);
isSP := function(M)       // M = +- permutation matrix ?
  for e in [F3!1, F3!-1] do
    N := e*M;
    if &and[ #[x : x in Eltseq(N[i]) | x ne 0] eq 1 and &+Eltseq(N[i]) eq 1 : i in [1..5] ] and
       &and[ #[i : i in [1..5] | N[i][j] ne 0] eq 1 : j in [1..5] ] then return true; end if;
  end for;
  return false;
end function;
GL5 := GL(5, F3);
Gaut := sub< GL5 | [GL5!M : M in auts] >;
Report(#auts eq 240 and #Gaut eq 240 and &and[ isSP(M) : M in auts ],
      "240 monomial automorphisms = (permutations of the 5 coordinates) x {+1,-1}");
idx := func< w | Index(Tcols, w) gt 0 select Index(Tcols, w) else Index(Tcols, -w) >;
perms := [ Sym(15)![ idx(Tcols[l]*M) : l in [1..15] ] : M in auts ];
G15 := sub< Sym(15) | perms >;
Report(IsIsomorphic(G15, Sym(5)) and Sort([#O : O in Orbits(G15)]) eq [5, 10] and
      { Multiplet(Cols[Rep(O)]) : O in Orbits(G15) | #O eq 10 } eq {"10"},
      "action on qutrits is S5 with orbits 10 and 5bar; full group S5 x Z2");
print "";
