"""Spherical coupled profiles (J-01): analytic coframe jets and Cartan curvature.

All coefficients are assembled before quadrature, never fitted in frequency.
Curvature is R0+n*R1+n**2*R2 at the base phase; n is a positive integer for
a smooth period-one profile. r,L,eta,V and nonzero amplitude are fixed.
"""
import itertools
import numpy as np

I=np.eye(4)
E3=np.zeros((3,3,3));E4=np.zeros((4,4,4,4))
for size,E in [(3,E3),(4,E4)]:
    for p in itertools.permutations(range(size)):
        E[p]=(-1)**sum(p[i]>p[j] for i in range(size) for j in range(i+1,size))
BASE=np.zeros((4,4,4));BASE[:3,:3,:3]=4*E3

def lc(C):return (C+C.transpose(2,1,0)-C.transpose(1,0,2))/2
def torsion(h,V):
    T=np.zeros((4,4,4));T[:3,:3,:3]=2*h*E3
    return T+V*(np.einsum('ac,b->abc',I,I[:,3])-np.einsum('ab,c->abc',I,I[:,3]))/3
def contortion(T):return -(T+T.transpose(2,0,1)-T.transpose(1,2,0))/2
def ricci(R):return np.einsum('abad->bd',R)
def tf(R):
    Q=ricci(R);sc=np.trace(Q)
    return R-(np.einsum('ac,bd->abcd',I,Q)-np.einsum('ad,bc->abcd',I,Q)-np.einsum('bc,ad->abcd',I,Q)+np.einsum('bd,ac->abcd',I,Q))/2+sc*(np.einsum('ac,bd->abcd',I,I)-np.einsum('ad,bc->abcd',I,I))/6
def project(R):
    S=(R+R.transpose(2,3,0,1))/2;F=(S+S.transpose(0,3,1,2)+S.transpose(0,2,3,1))/3
    W=tf(S-F);TF=tf(R)
    return {'W':W,'R2':TF-W-F,'R3':F,'TF':TF}
def product(G,H):return np.einsum('aec,ebd->abcd',G,H)-np.einsum('aed,ebc->abcd',G,H)
def structure_term(G,C):return np.einsum('abe,ecd->abcd',G,C)
def derivative_term(Gp,B):return np.einsum('c,abd->abcd',B[3],Gp)-np.einsum('d,abc->abcd',B[3],Gp)

def coframe_jets(y,amplitude=.5,r=3.,L=1.,frequency=1.):
    if not (r>0 and L>0 and 0<=amplitude<np.sqrt(2)*r/L):raise ValueError('nondegenerate fixed-amplitude X4 family required')
    phase=2*np.pi*frequency*y
    w=amplitude*np.cos(phase);wp=-2*np.pi*frequency*amplitude*np.sin(phase);wpp=-(2*np.pi*frequency)**2*w
    b=L*L/(2*r*r);F=1-b*w*w;s=.5*np.log(F)
    sp=-b*w*wp/F;spp=-b*((1+b*w*w)*wp*wp/F**2+w*wpp/F)
    A=np.diag([r*np.exp(s),r*np.exp(-s),r,L]);A[3,0]=L*w
    Ap=np.diag([A[0,0]*sp,-A[1,1]*sp,0.,0.]);Ap[3,0]=L*wp
    App=np.diag([A[0,0]*(sp*sp+spp),A[1,1]*(sp*sp-spp),0.,0.]);App[3,0]=L*wpp
    return A,Ap,App,{'w':w,'wp':wp,'s':s,'sp':sp,'M':np.diag([np.exp(s),np.exp(-s),1.])}

def structures(A,Ap,App):
    B=np.linalg.inv(A);Bp=-B@Ap@B
    C0=np.einsum('ai,jb,kc,ijk->abc',A,B,B,BASE)
    C0p=(np.einsum('ai,jb,kc,ijk->abc',Ap,B,B,BASE)+np.einsum('ai,jb,kc,ijk->abc',A,Bp,B,BASE)+np.einsum('ai,jb,kc,ijk->abc',A,B,Bp,BASE))
    Q=np.einsum('af,b,fc->abc',Ap,B[3],B)
    Qp=np.einsum('af,b,fc->abc',App,B[3],B)+np.einsum('af,b,fc->abc',Ap,Bp[3],B)+np.einsum('af,b,fc->abc',Ap,B[3],Bp)
    return C0,Q-Q.transpose(0,2,1),C0p,Qp-Qp.transpose(0,2,1),B

def frequency_curvatures(y,amplitude=.5,r=3.,L=1.,eta=1.,V=1.):
    """Return R(n)=R0+n R1+n^2 R2 with coefficient functions at phase 2 pi y.

    For integer n the complete circle repeats these coefficients n times.
    Derivatives of contortion vanish in the chosen frame since r,eta,V are fixed;
    its connection and noncoordinate products are retained in every coefficient.
    """
    A,Ap,App,info=coframe_jets(y,amplitude,r,L)
    C0,C1,C0p,C1p,B=structures(A,Ap,App)
    LC0,LC1=lc(C0),lc(C1);G0p,G1p=lc(C0p),lc(C1p)
    K=contortion(torsion(eta/r,V))
    def coefficients(G0):
        return np.array([
            product(G0,G0)+structure_term(G0,C0),
            product(G0,LC1)+product(LC1,G0)+structure_term(G0,C1)+structure_term(LC1,C0)+derivative_term(G0p,B),
            product(LC1,LC1)+structure_term(LC1,C1)+derivative_term(G1p,B)])
    return coefficients(LC0+K),coefficients(LC0),info,(C0,C1,LC0+K,LC1)

def quadratic_coefficients(tensors):
    return np.array([sum(np.sum(tensors[i]*tensors[j]) for i in range(3) for j in range(3) if i+j==degree) for degree in range(5)])

def invariants_at_phase(y,**kwargs):
    R,RL,info,aux=frequency_curvatures(y,**kwargs)
    P=[project(x) for x in R];PL=[project(x) for x in RL]
    vals={name:quadratic_coefficients(np.array([p[key] for p in P])) for name,key in [('EC_IRR','W'),('EC_TRACEFREE','TF'),('R2','R2'),('R3','R3')]}
    vals['LC']=quadratic_coefficients(np.array([p['W'] for p in PL]))
    vals['scalar']=np.array([np.trace(ricci(x)) for x in R])
    vals['scalar_LC']=np.array([np.trace(ricci(x)) for x in RL])
    return vals,info

def integrate(samples=512,**kwargs):
    """Average analytic frequency coefficients over an equally spaced circle grid."""
    if not isinstance(samples,(int,np.integer)) or samples<1:raise ValueError('samples must be a positive integer')
    values=[invariants_at_phase(j/samples,**kwargs)[0] for j in range(samples)]
    return {key:np.mean([v[key] for v in values],axis=0).tolist() for key in values[0]}

def action_polynomial(mean,*,model='EC_IRR',alpha=0.,kappa=1.,r=3.,L=1.):
    """Return the five ascending coefficients of the integrated circle action.

    mean is returned by integrate with the same r,L. Complete NY is zero.
    Odd averages vanish exactly by reflection; quadrature roundoff is removed.
    """
    if model not in ('LC','EC_IRR','EC_TRACEFREE'):raise ValueError('Unknown curvature model')
    if kappa<=0 or r<=0 or L<=0:raise ValueError('kappa,r,L must be positive')
    volume=np.pi**2*L*r**3/4
    coefficients=-alpha*volume*np.asarray(mean[model])
    coefficients[:3]-=volume*np.asarray(mean['scalar'])/(2*kappa*kappa)
    coefficients[[1,3]]=0
    return coefficients
