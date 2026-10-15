"""Monomial automorphism group of the fifteen-qudit code for p = 5, 7, 11:
S_5 x F_p^*, of order 120 (p - 1).

Counts projective linear maps that permute the 15 columns up to scalars and preserve the
stabilizer space; the monomial group has order (p-1) times this count. Also checks that the 120
permutations of the five coordinates of q are among them, so they are all of them. Uses the generator matrix
printed in the paper. Asserts a count of 120 for every p tested. Run time about one minute per prime.
The p = 3 case is automorphisms_15.py; Magma script 06 checks all four primes independently.
"""
import itertools, numpy as np, sys, time
L=[-1,-1,-1,-1,1,1,1,1,1,1, 1,-1,-1,-1,-1]
S=[[0,-1,-1,-1,1,1,1,0,0,0, 1,-1,0,0,0],[-1,1,0,0,0,-1,-1,1,1,0, 0,1,-1,0,0],
   [0,-1,1,0,-1,1,0,0,-1,1, 0,0,1,-1,0],[0,0,-1,1,0,-1,1,-1,1,0, 0,0,0,1,-1]]
def inv(M,p):
    n=len(M); A=[list(map(int,r))+[int(i==j) for j in range(n)] for i,r in enumerate(M)]
    for c in range(n):
        r=next((i for i in range(c,n) if A[i][c]%p),None)
        if r is None: return None
        A[c],A[r]=A[r],A[c]; iv=pow(A[c][c],-1,p); A[c]=[x*iv%p for x in A[c]]
        for i in range(n):
            if i!=c and A[i][c]%p: f=A[i][c]; A[i]=[(a-f*b)%p for a,b in zip(A[i],A[c])]
    return np.array([r[n:] for r in A])
def run(p):
    H=np.array([L]+S)%p; cols=[H[:,l] for l in range(15)]
    def norm(v):
        v=v%p; i=next(i for i in range(5) if v[i]); return tuple(v*pow(int(v[i]),-1,p)%p)
    pts={norm(c):l for l,c in enumerate(cols)}
    # frame: basis + a point with all-nonzero coordinates in that basis
    for basis in itertools.combinations(range(15),5):
        F=np.array([cols[l] for l in basis]).T
        Fi=inv(F,p)
        if Fi is None: continue
        for l6 in range(15):
            if l6 in basis: continue
            c=Fi.dot(cols[l6])%p
            if all(c): frame=list(basis)+[l6]; coef=c; break
        else: continue
        break
    Fi=inv(np.array([cols[l] for l in frame[:5]]).T,p)
    count=0
    for img in itertools.permutations(range(15),5):
        G=np.array([cols[l] for l in img]).T
        Gi=inv(G,p)
        if Gi is None: continue
        for l6 in range(15):
            if l6 in img: continue
            d=Gi.dot(cols[l6])%p
            if not all(d): continue
            # B f_i = mu_i g_i with mu_i = d_i/coef_i (up to global scalar)
            mu=[d[i]*pow(int(coef[i]),-1,p)%p for i in range(5)]
            B=(G.dot(np.diag(mu)).dot(Fi))%p
            imgs=[pts.get(norm(B.dot(c)%p)) for c in cols]
            if None in imgs or len(set(imgs))<15: continue
            # H M = B H: stabilizer rows map to stabilizer rows iff B[1:, 0] = 0
            if any(B[1:,0]%p): continue
            count+=1
    return count
def coordinate_permutations_preserve(p):
    """The 120 permutations of the five coordinates of q preserve the row space and the stabilizer
    space, so with the count above they are all of the projective automorphisms."""
    from sm16 import rank_mod
    W = [tuple(-1 if k in ij else 1 for k in range(5)) for ij in itertools.combinations(range(5), 2)]
    W += [tuple(1 if t == k else -1 for t in range(5)) for k in range(5)]
    H = np.array([L] + S) % p
    for P in itertools.permutations(range(5)):
        perm = [W.index(tuple(q[P[i]] for i in range(5))) for q in W]
        M = H[:, perm]
        if rank_mod(np.vstack([H, M]).tolist(), p) != 5 or rank_mod(np.vstack([H[1:], M[1:]]).tolist(), p) != 4:
            return False
    return True
if __name__ == '__main__':
    for p in [5,7,11]:
        t=time.time(); n=run(p)
        print('p=%d: projective automorphisms %d, monomial group order %d = (p-1)*%d'%(p,n,n*(p-1),n),
              ' time %.0fs'%(time.time()-t)); sys.stdout.flush()
        assert n == 120, (p, n)
        assert coordinate_permutations_preserve(p)
    print('monomial automorphism group S_5 x F_p^* for p = 5, 7, 11')
