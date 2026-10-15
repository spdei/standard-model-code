// common.m -- shared definitions for the verification scripts.
// Conventions: a fermion is its charge vector q in {+1,-1}^5 (weights of the 16 of so(10),
// scaled by 2).  The 15 Standard Model fields are listed in the column order of the
// The paper's table of the sixteen weights as Standard Model fields (the doubling order), so all
// scripts use the same column order.

Z := Integers(); Q := Rationals();

// Weights of a spinor of so(2m): sign vectors in {+-1}^m with an even number of minus signs.
function SpinorWeights(m)
  W := [];
  for i in [0..2^m-1] do
    b := Intseq(i, 2); b cat:= [0 : k in [1..m-#b]];
    q := [1 - 2*x : x in b];
    if #[x : x in q | x eq -1] mod 2 eq 0 then Append(~W, q); end if;
  end for;
  return W;
end function;

Nu := [1,1,1,1,1];                       // nu^c, the removed field

// Column order of the paper's table of fields (= column order of H15), with field names.
Names := ["dc1","dc2","uc3","dc3","uc2","uc1","ec",
          "d1","d2","u3","d3","u2","u1","e","nu"];
Cols  := [[ 1,-1,-1,-1,-1], [-1, 1,-1,-1,-1], [-1,-1, 1, 1, 1], [-1,-1, 1,-1,-1],
          [-1, 1,-1, 1, 1], [ 1,-1,-1, 1, 1], [ 1, 1, 1,-1,-1],
          [-1, 1, 1,-1, 1], [ 1,-1, 1,-1, 1], [ 1, 1,-1, 1,-1], [ 1, 1,-1,-1, 1],
          [ 1,-1, 1, 1,-1], [-1, 1, 1, 1,-1], [-1,-1,-1,-1, 1], [-1,-1,-1, 1,-1]];
All16 := Cols cat [Nu];                  // index 16 is nu^c

// Multiplet of a field: "10" (two minus signs) or "5bar" (four minus signs).
Multiplet := func< q | #[x : x in q | x eq -1] eq 2 select "10" else "5bar" >;

// Componentwise product of vectors.
CP := func< u, v | Vector([u[i]*v[i] : i in [1..Degree(u)]]) >;

// Integer trilinear form with weights c.
TauZ := func< x, y, z, c | &+[Z| c[l]*x[l]*y[l]*z[l] : l in [1..#c]] >;

// Equation (dist): smallest |J| with s in span{V[j] : j in J}, searching |J| <= maxw.
// Returns 0 if no such J with |J| <= maxw exists.
function ZDistVec(V, s, maxw)
  n := #V;
  for w in [1..maxw] do
    for J in Subsets({1..n}, w) do
      M := Matrix([V[j] : j in J]);
      if Rank(M) eq Rank(VerticalJoin(M, Matrix([s]))) then return w; end if;
    end for;
  end for;
  return 0;
end function;

// Same, for integer columns reduced modulo p.
ZDistModP := func< cols, s, p, maxw |
  ZDistVec([Vector(GF(p), c) : c in cols], Vector(GF(p), s), maxw) >;

procedure Report(cond, msg)
  if cond then printf "  [ok]   %o\n", msg;
  else printf "  [FAIL] %o\n", msg; end if;
end procedure;
