"""Exact derivations for the restricted families in Sections 3--5."""
import itertools
import sympy as s
from .cartan import *
from . import exterior as f

def derive():
    checks={}
    def ck(name,value):checks[name]=str(clean(value))
    def cn(name,A):ck(name,sum(v*v for v in s.flatten(A.tolist())))
    n=s.symbols('n0:3',real=True)
    h,V,dh,dV,ddh,ddV=s.symbols('h V h_dot V_dot h_ddot V_ddot',real=True)
    def D(i,x):return (s.diff(x,h)*dh+s.diff(x,V)*dV+s.diff(x,dh)*ddh+s.diff(x,dV)*ddV) if i==3 else Z
    C=structure(n);T=torsion(h,V);LC,G=connection(C,T)
    R=curvature(G,C,D);RL=curvature(LC,C);PP=project(R);PL=project(RL)
    vals={k:norm(v) for k,v in PP.items()};vals.update(densities(R,T))
    RLC=clean(s.trace(ricci(RL)));I0=norm(PL['W'])
    ck('circle_W_tensor',norm(arr(4,lambda *i:PP['W'][i]-PL['W'][i])))
    ck('circle_R2',vals['R2']-6*(dh+2*h*V/3)**2)
    ck('circle_R3',vals['R3']-6*(dh-2*h*V/3)**2)
    ck('circle_TF',vals['TF']-I0-12*dh**2-s.Rational(16,3)*h*h*V*V)
    ck('circle_scalar',vals['R']-RLC+6*h*h+s.Rational(2,3)*V*V-2*dV)
    ck('circle_NY',vals['NY']+6*dh)
    ck('circle_P_transgression',vals['P']-D(3,2*RLC*h-4*h**3-s.Rational(4,3)*h*V*V))
    de=f.de_from_C(C);om=f.connection_forms(G);F=f.curvature_forms(om,de,D);RF=f.R_from_forms(F)
    ck('curvature_vs_exterior_forms',norm(arr(4,lambda *i:R[i]-RF[i])))
    ck('P_vs_wedge',vals['P']-f.metrics(R)['wedge_external'])
    ck('Chern_Simons',f.exterior(f.cs_trace(om,de,D),de,D).get(f.TOP,0)+vals['P'])
    e=[{(a,):s.S.One} for a in range(4)]
    Tf=[{(b,c):T[a,b,c] for b,c in f.PAIRS if T[a,b,c]!=0} for a in range(4)]
    for a in range(4):
        first=f.add(de[a],*[f.wedge(om[a][b],e[b]) for b in range(4)],f.scale(-1,Tf[a]))
        ck('first_Cartan_'+str(a),sum(v*v for v in first.values()))
        b1=f.add(f.exterior(Tf[a],de,D),*[f.wedge(om[a][b],Tf[b]) for b in range(4)],*[f.scale(-1,f.wedge(F[a][b],e[b])) for b in range(4)])
        ck('first_Bianchi_'+str(a),sum(v*v for v in b1.values()))
        for b in range(4):
            b2=f.add(f.exterior(F[a][b],de,D),*[f.wedge(om[a][c],F[c][b]) for c in range(4)],*[f.scale(-1,f.wedge(F[a][c],om[c][b])) for c in range(4)])
            ck('second_Bianchi_'+str(a)+str(b),sum(v*v for v in b2.values()))
    ny=f.exterior(f.add(*[f.wedge(e[a],Tf[a]) for a in range(4)]),de,D).get(f.TOP,0)
    ck('NY_vs_exterior',ny-vals['NY'])
    for a,b in [('W','R2'),('W','R3'),('R2','R3')]:ck(a+'_'+b+'_orthogonal',inner(PP[a],PP[b]))
    pair=list(itertools.combinations(range(4),2));maps={k:[] for k in PP}
    for ii,jj in itertools.product(range(6),repeat=2):
        basis=s.MutableDenseNDimArray.zeros(4,4,4,4);a,b=pair[ii];c,d=pair[jj]
        basis[a,b,c,d]=basis[b,a,d,c]=s.S.One
        basis[b,a,c,d]=basis[a,b,d,c]=-s.S.One
        bp=project(basis)
        for key in maps:maps[key].append([bp[key][*p,*q] for p in pair for q in pair])
    maps={k:s.Matrix(36,36,lambda i,j:columns[j][i]) for k,columns in maps.items()}
    for key,rank in [('W',10),('R2',9),('R3',1),('TF',20)]:
        mm=maps[key];ck(key+'_rank',mm.rank()-rank);cn(key+'_idempotent',mm*mm-mm);cn(key+'_selfadjoint',mm-mm.T)
    cn('projector_reconstruction',maps['TF']-maps['W']-maps['R2']-maps['R3'])
    for a,b in [('W','R2'),('W','R3'),('R2','R3')]:cn(a+'_'+b+'_orthogonal_maps',maps[a]*maps[b])
    lam=[(n[i]**2-(n[(i+1)%3]-n[(i+2)%3])**2)/2 for i in range(3)]
    ck('LC_Ricci_eigenvalues',RLC-sum(lam));ck('LC_Weyl_normalization',I0-2*sum(v*v for v in lam)+2*RLC*RLC/3)
    ck('projected_entry_R2',PP['R2'][0,1,2,3]-dh/2-h*V/3)
    ck('projected_entry_R3',PP['R3'][0,1,2,3]-dh/2+h*V/3)
    # Pure axial first-jet control in all four orthonormal directions.
    g=s.symbols('g0:4',real=True)
    _,GA=connection(C,torsion(h,0));RA=curvature(GA,C,lambda i,x:s.diff(x,h)*g[i]);PA=project(RA)
    ck('AX_spatial_first_jet',norm(PA['TF'])-I0-4*sum(t*t for t in g[:3])-12*g[3]**2)
    # Constant twist e3 -> e3 + w_i e_i, not a varying twist profile.
    w=s.symbols('w0:3',real=True);CT=structure(n)
    for j,k in itertools.product(range(3),repeat=2):CT[3,j,k]=sum(w[i]*CT[i,j,k] for i in range(3))
    LT,GT=connection(CT,T);RT=curvature(GT,CT);PT=project(RT);WT=project(curvature(LT,CT))['W']
    twist={k:norm(v) for k,v in PT.items()};twist.update(densities(RT,T))
    ck('twist_W_tensor',norm(arr(4,lambda *i:PT['W'][i]-WT[i])))
    ck('twist_R2',twist['R2']-s.Rational(8,3)*h*h*V*V)
    ck('twist_R3',twist['R3']-s.Rational(8,3)*h*h*V*V)
    ck('twist_NY',twist['NY']);ck('twist_P',twist['P'])
    # Torsion axis u in the fixed product frame. Set |u|=1 for a rotation.
    u=s.symbols('u0:4',real=True);u2=sum(t*t for t in u)
    TU=arr(3,lambda a,b,c:2*h*sum(s.LeviCivita(a,b,c,d)*u[d] for d in range(4))+V*(delta(a,c)*u[b]-delta(a,b)*u[c])/3)
    _,GU=connection(C,TU);RU=curvature(GU,C);PU=project(RU)
    orient={k:norm(v) for k,v in PU.items()};orient.update(densities(RU,TU))
    E=sum((n[(i+1)%3]-n[(i+2)%3])**2*u[i]**2 for i in range(3))
    ck('orientation_W_tensor',norm(arr(4,lambda *i:PU['W'][i]-PL['W'][i])))
    ck('orientation_R2',orient['R2']-s.Rational(8,3)*h*h*V*V*u2*u2-4*h*h*E)
    ck('orientation_R3',orient['R3']-s.Rational(8,3)*h*h*V*V*u2*u2)
    ck('orientation_scalar',orient['R']-RLC+(6*h*h+s.Rational(2,3)*V*V)*u2)
    ck('orientation_NY',orient['NY']);ck('orientation_P',orient['P'])
    # Exact Hessians in the normalized exponential shape chart.
    t,r,L,kappa,alpha=s.symbols('t r L kappa alpha',positive=True)
    CS=structure([4*t*t/r,4*t*t/r,4/(t**4*r)]);LS,_=connection(CS,torsion(0,0));RS=curvature(LS,CS)
    DS=lambda x:s.diff(x,t)*t/s.sqrt(6)
    ck('S3_R_shape_second',DS(DS(s.trace(ricci(RS)))).subs(t,1)+64/r**2)
    ck('S3_W_shape_second',DS(DS(norm(project(RS)['W']))).subs(t,1)-4096/r**4)
    A,B=s.symbols('A B',real=True);fp=2*A+4*B
    HN=s.diag((32*A+128*B)/3,0,0,4*A+8*B,4*A+8*B)
    shape_basis=[s.diag(1,1,-2)/s.sqrt(6),s.diag(1,-1,0)/s.sqrt(2),s.Matrix([[0,1,0],[1,0,0],[0,0,0]])/s.sqrt(2),s.Matrix([[0,0,1],[0,0,0],[1,0,0]])/s.sqrt(2),s.Matrix([[0,0,0],[0,0,1],[0,1,0]])/s.sqrt(2)]
    za=s.Matrix([2*t[2,2] for t in shape_basis])
    zab=s.Matrix(5,5,lambda i,j:2*(shape_basis[i]*shape_basis[j]+shape_basis[j]*shape_basis[i])[2,2])
    cn('Nil_full_shape_H',(2*A+12*B)*za*za.T+fp*zab-HN)
    eps,sh=s.symbols('epsilon s',real=True);z=(1+eps)**s.Rational(-4,3)/(1+sh)**2
    oldH=s.hessian(A*z*z+B*z**4,[eps,sh]).subs({eps:0,sh:0})
    ck('Nil_old_cross',oldH[0,1]-s.Rational(32,3)*(A+4*B))
    logf=s.Rational(2,3)*s.log(1+eps)+s.log(1+sh)
    q=s.Matrix([s.sqrt(6)*logf/2,-logf/s.sqrt(2),0,0,0]);xx=[eps,sh];origin={eps:0,sh:0}
    J=q.jacobian(xx).subs(origin);grad=fp*za
    mapped=J.T*HN*J+s.Matrix(2,2,lambda a,b:sum(grad[i]*s.diff(q[i],xx[a],xx[b]).subs(origin) for i in range(5)))
    ck('Nil_old_tangent_rank',J.rank()-1);cn('Nil_coordinate_H_with_gradient',mapped-oldH)
    # Euler equations: alpha here is arbitrary real; positivity is not used.
    aa,kk=s.symbols('alpha kappa',real=True,nonzero=True)
    ell=3*h*h/kk**2+V*V/(3*kk**2)-aa*(12*dh*dh+s.Rational(16,3)*h*h*V*V)
    euler_h=clean(s.diff(ell,h)-D(3,s.diff(ell,dh)));euler_V=clean(s.diff(ell,V))
    ck('Euler_h',euler_h-(6/kk**2-s.Rational(32,3)*aa*V*V)*h-24*aa*ddh)
    ck('Euler_V',euler_V-(2/(3*kk**2)-s.Rational(32,3)*aa*h*h)*V)
    # Static sections of the general circle result on all three geometries.
    stat={}
    for name,nn in [('S3',[4/r]*3),('T3',[0]*3),('Nil3',[0,0,-1/r])]:
        sub=dict(zip(n,nn));sub.update({dh:0,dV:0})
        stat[name]={k:str(clean(v.subs(sub))) for k,v in vals.items()}
    failed={k:v for k,v in checks.items() if v!='0'}
    if failed:raise ArithmeticError(f'Exact identity residuals: {failed}')
    return {'schema':'paper03ec.exact.v1','conventions':'Euclidean; T geometric; circle index 3; h=eta/r; dot=d/dz; y period 1; z=L*y',
            'circle':{k:str(v) for k,v in vals.items()},'R_LC':str(RLC),'I_LC':str(I0),
            'twist':{k:str(v) for k,v in twist.items()},'orientation':{k:str(v) for k,v in orient.items()},
            'static':stat,'Nil_gradient':[str(-4*fp/s.sqrt(6)),'0','0','0','0'],
            'Nil_H':[[str(v) for v in HN.row(i)] for i in range(5)],'Nil_old_path_H':[[str(v) for v in oldH.row(i)] for i in range(2)],
            'Euler_h':str(euler_h),'Euler_V':str(euler_V),'identity_residuals':checks}
