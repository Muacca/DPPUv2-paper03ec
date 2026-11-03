"""Exterior forms for independently reproducing the paper's Cartan identities."""
from itertools import combinations, permutations, product
import sympy as s
PAIRS=list(combinations(range(4),2))
TOP=(0,1,2,3)
def clean(f):
    return {k:v for k,w in f.items() if (v:=s.cancel(w))!=0}
def add(*fs):
    out={}
    for f in fs:
        for k,v in f.items():out[k]=out.get(k,0)+v
    return clean(out)
def scale(c,f):return clean({k:c*v for k,v in f.items()})
def wedge(f,g):
    out={}
    for a,v in f.items():
        for b,w in g.items():
            seq=a+b
            if len(set(seq))<len(seq):continue
            sign=(-1)**sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))
            key=tuple(sorted(seq))
            out[key]=out.get(key,0)+sign*v*w
    return clean(out)
def exterior(f,de,deriv=None):
    out={}
    for inds,c in f.items():
        if deriv is not None:
            out=add(out,wedge({(i,):deriv(i,c) for i in range(4)},{inds:s.S.One}))
        for j,a in enumerate(inds):
            out=add(out,scale((-1)**j*c,wedge(wedge({inds[:j]:s.S.One},de[a]),{inds[j+1:]:s.S.One})))
    return clean(out)
def de_from_C(C):return [{(b,c):C[a,b,c] for b,c in PAIRS if C[a,b,c]!=0} for a in range(4)]
def connection_forms(G):return [[{(c,):G[a,b,c] for c in range(4) if G[a,b,c]!=0} for b in range(4)] for a in range(4)]
def curvature_forms(omega,de,deriv=None):
    return [[add(exterior(omega[a][b],de,deriv),*[wedge(omega[a][c],omega[c][b]) for c in range(4)]) for b in range(4)] for a in range(4)]
def R_from_forms(F):
    R=s.MutableDenseNDimArray.zeros(4,4,4,4)
    for a,b in product(range(4),repeat=2):
        for (c,d),val in F[a][b].items():R[a,b,c,d]=val;R[a,b,d,c]=-val
    return R
def forms_from_R(R):return [[{(c,d):R[a,b,c,d] for c,d in PAIRS if R[a,b,c,d]!=0} for b in range(4)] for a in range(4)]
def pfaff(v):return s.cancel(v[0]*v[5]-v[1]*v[4]+v[2]*v[3])
def metrics(R):
    M=s.Matrix(6,6,lambda i,j:R[*PAIRS[i],*PAIRS[j]])
    internal=s.cancel(sum(pfaff(list(M[:,j])) for j in range(6)))
    external=s.cancel(4*sum(pfaff(list(M[i,:])) for i in range(6)))
    q_internal=s.cancel(8*internal)
    q_external=s.cancel(2*external)
    F=forms_from_R(R)
    true_wedge=add(*[wedge(F[a][b],F[a][b]) for a,b in product(range(4),repeat=2)]).get(TOP,0)
    euler=add(*[scale(s.LeviCivita(a,b,c,d),wedge(F[a][b],F[c][d])) for a,b,c,d in permutations(range(4))]).get(TOP,0)
    return {'block_internal':internal,'pontryagin_external':external,'wedge_external':true_wedge,
            'full_internal_inner':q_internal,'full_external_inner':q_external,'euler_numerator':euler,
            'column_pfaffians':[pfaff(list(M[:,j])) for j in range(6)],
            'row_pfaffians':[pfaff(list(M[i,:])) for i in range(6)],'pair_matrix':M}
def cs_trace(omega,de,deriv=None):
    # Matrix trace CS, so dCS=Tr(F wedge F)=-sum_ab F_ab wedge F_ab.
    quadratic=add(*[wedge(omega[a][b],exterior(omega[b][a],de,deriv)) for a,b in product(range(4),repeat=2)])
    cubic=add(*[wedge(wedge(omega[a][b],omega[b][c]),omega[c][a]) for a,b,c in product(range(4),repeat=3)])
    return add(quadratic,scale(s.Rational(2,3),cubic))
def serialized(obj):
    if isinstance(obj,dict):return {str(k):serialized(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):return [serialized(v) for v in obj]
    if isinstance(obj,s.MatrixBase):return serialized(obj.tolist())
    return str(s.factor(obj))
