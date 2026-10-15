// 06_builtin_automorphisms.m -- monomial automorphism groups via Magma's built-in (Leon)
// automorphism groups of linear codes, as stated in
// the paper. The monomial automorphisms of the code are those of the row space H that also
// preserve the stabilizer space A, i.e. Aut(H) meet Aut(A).
//   p = 3:          order 240 = |S5 x Z2|, an independent cross-check of script 05;
//   p = 5, 7, 11:   order 120 (p-1) = |S5 x F_p^*|, the Weyl group fixing nu^c and global scalars.
load "common.m";
print "=== 06: automorphism groups via built-in AutomorphismGroup ===";
for p in [3, 5, 7, 11] do
  F := GF(p);
  H := LinearCode(Matrix(F, [[q[i] : q in Cols] : i in [1..5]]));
  A := LinearCode(Matrix(F, [[q[i] - q[i+1] : q in Cols] : i in [1..4]]));
  Report(Dimension(H) eq 5 and Dimension(A) eq 4 and A subset H,
         Sprintf("p = %o: H = [15,5]_%o row space, A = [15,4]_%o stabilizers", p, p, p));
  expected := 120 * (p - 1);
  try
    AutH := AutomorphismGroup(H);
    AutA := AutomorphismGroup(A);
    G := AutH meet AutA;
    printf "    |Aut(H)| = %o, |Aut(A)| = %o, |Aut(H) meet Aut(A)| = %o\n", #AutH, #AutA, #G;
    Report(#G eq expected,
           Sprintf("p = %o: Aut(H) meet Aut(A) has order %o = |S5 x F_%o^*|", p, expected, p));
  catch err
    printf "  [FAIL] p = %o: built-in automorphism computation raised an error:\n", p;
    print err;
  end try;
end for;
print "";
