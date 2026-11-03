"""Four-dimensional Euclidean Cartan tensors used in the paper.

Indices 0,1,2 are spatial and 3 is the circle. T is geometric torsion.
Curvature includes derivatives only when the explicit deriv callback is supplied.
The orthogonal projections are W (10), R2 (9), R3 (1), and their sum TF (20).
"""
import itertools
import sympy as s
IX3=list(itertools.product(range(4),repeat=3))
IX4=list(itertools.product(range(4),repeat=4))
Z=s.S.Zero

def arr(rank,fn):
    return s.MutableDenseNDimArray([fn(*i) for i in itertools.product(range(4),repeat=rank)],(4,)*rank)

def clean(x): return s.factor(s.cancel(x))
def delta(a,b): return s.S.One if a==b else Z

def structure(n):
    return arr(3,lambda a,b,c:n[a]*s.LeviCivita(a,b,c) if max(a,b,c)<3 else Z)

def transformed(C,A):
    inv=A.inv()
    return arr(3,lambda a,b,c:clean(sum(A[a,i]*inv[j,b]*inv[k,c]*C[i,j,k] for i,j,k in IX3 if A[a,i]!=0 and inv[j,b]!=0 and inv[k,c]!=0 and C[i,j,k]!=0)))

def torsion(h,V):
    return arr(3,lambda a,b,c:(2*h*s.LeviCivita(a,b,c) if max(a,b,c)<3 else Z)+V*(delta(a,c)*delta(b,3)-delta(a,b)*delta(c,3))/3)

def connection(C,T):
    LC=arr(3,lambda a,b,c:(C[a,b,c]+C[c,b,a]-C[b,a,c])/2)
    K=arr(3,lambda a,b,c:-(T[a,b,c]+T[b,c,a]-T[c,a,b])/2)
    return LC,arr(3,lambda a,b,c:LC[a,b,c]+K[a,b,c])

def curvature(G,C,deriv=None):
    def component(a,b,c,d):
        value=sum(G[a,e,c]*G[e,b,d]-G[a,e,d]*G[e,b,c]+G[a,b,e]*C[e,c,d] for e in range(4))
        if deriv is not None: value+=deriv(c,G[a,b,d])-deriv(d,G[a,b,c])
        return clean(value)
    return arr(4,component)

def ricci(R): return s.Matrix(4,4,lambda b,d:sum(R[a,b,a,d] for a in range(4)))

def tracefree(R):
    Q=ricci(R); sc=s.trace(Q)
    return arr(4,lambda a,b,c,d:clean(R[a,b,c,d]-(delta(a,c)*Q[b,d]-delta(a,d)*Q[b,c]-delta(b,c)*Q[a,d]+delta(b,d)*Q[a,c])/2+sc*(delta(a,c)*delta(b,d)-delta(a,d)*delta(b,c))/6))

def project(R):
    R=arr(4,lambda *i:s.sympify(R[i]))
    S=arr(4,lambda a,b,c,d:(R[a,b,c,d]+R[c,d,a,b])/2)
    F=arr(4,lambda a,b,c,d:clean((S[a,b,c,d]+S[a,c,d,b]+S[a,d,b,c])/3))
    W=tracefree(arr(4,lambda a,b,c,d:S[a,b,c,d]-F[a,b,c,d]))
    TF=tracefree(R)
    A=arr(4,lambda a,b,c,d:clean(TF[a,b,c,d]-W[a,b,c,d]-F[a,b,c,d]))
    return {'W':W,'R2':A,'R3':F,'TF':TF}

def inner(A,B): return clean(sum(A[i]*B[i] for i in IX4))
def norm(A): return inner(A,A)

def densities(R,T):
    eps=[(p,s.LeviCivita(*p)) for p in itertools.permutations(range(4))]
    tt=clean(sum(e*sum(T[i,a,b]*T[i,c,d] for i in range(4))/4 for (a,b,c,d),e in eps))
    ree=clean(sum(e*R[a,b,c,d]/2 for (a,b,c,d),e in eps))
    p=clean(sum(e*sum(R[a,b,c,d]*R[a,b,f,g] for a in range(4) for b in range(4))/4 for (c,d,f,g),e in eps))
    return {'R':clean(s.trace(ricci(R))),'TT':tt,'REE':ree,'NY':clean(tt-ree),'P':p}

def static(topology,r,L,M=None,w=None,O=None,h=Z,V=Z):
    M=s.eye(3) if M is None else M
    w=s.zeros(3,1) if w is None else s.Matrix(w)
    O=s.eye(4) if O is None else O
    A=s.diag(r,r,r,L);A[:3,:3]=r*M;A[3,:3]=(L*w).T
    A=O*A
    C=transformed(structure({'S3':[4,4,4],'T3':[0,0,0],'Nil3':[0,0,-1]}[topology]),A)
    T=torsion(h,V);LC,G=connection(C,T)
    R=curvature(G,C);RLC=curvature(LC,C)
    p=project(R)
    return {'A':A,'C':C,'T':T,'LC':LC,'G':G,'R':R,'RLC':RLC,'projected':p,'norms':{k:norm(v) for k,v in p.items()},'densities':densities(R,T)}
