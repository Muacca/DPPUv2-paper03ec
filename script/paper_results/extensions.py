"""Exact stationary, coupled-profile, and fixed-quotient formulas.

The proofs of completeness, positivity, and the infinite-lattice minimum are
in Sections 3, 5, 6 and Appendix E. Here tensor algebra reproduces their
coefficients. No external validation directory or prior-paper package is used.
"""
import itertools
import sympy as s
from . import cartan as c

def stationary_formulas():
    checks={}
    def ck(name,value):checks[name]=str(c.clean(s.simplify(value)))
    k,alpha,L=s.symbols('kappa alpha L',positive=True)
    h,V,hd,hdd,phi,v,pd=s.symbols('h V h_dot h_ddot phi v phi_dot',real=True)
    ell=3*h*h/k**2+V*V/(3*k**2)-alpha*(12*hd*hd+s.Rational(16,3)*h*h*V*V)
    Eh=s.diff(ell,h)-s.diff(s.diff(ell,hd),hd)*hdd
    EV=s.diff(ell,V)
    H=1/(4*k*s.sqrt(alpha));K=3*H
    for sh,sv in itertools.product([-1,1],repeat=2):
        sub={h:sh*H,V:sv*K,hd:0,hdd:0}
        ck(f'mixed Euler_h {sh},{sv}',Eh.subs(sub));ck(f'mixed Euler_V {sh},{sv}',EV.subs(sub))
        ck(f'mixed action {sh},{sv}',ell.subs(sub)-3/(16*alpha*k**4))
    z=s.symbols('z',real=True);n=s.symbols('n',positive=True,integer=True);a,b=s.symbols('a b',real=True)
    an=L**2/(16*s.pi**2*k**2*n*n);omega=2*s.pi*n/L
    profile=a*s.cos(omega*z)+b*s.sin(omega*z)
    ck('resonant harmonic Euler',Eh.subs({alpha:an,h:profile,V:0,hdd:s.diff(profile,z,2)}))
    ck('resonant harmonic mean action',s.integrate(ell.subs({alpha:an,h:profile,V:0,hd:s.diff(profile,z)}),(z,0,L)))
    eps=s.symbols('eps',real=True)
    second=s.diff(ell.subs({h:h+eps*phi,V:V+eps*v,hd:hd+eps*pd}),eps,2).subs(eps,0)
    expected=(6/k**2-s.Rational(32,3)*alpha*V*V)*phi**2+(2/(3*k**2)-s.Rational(32,3)*alpha*h*h)*v*v-s.Rational(128,3)*alpha*h*V*phi*v-24*alpha*pd*pd
    ck('full restricted second variation',second-expected)
    xi=s.symbols('xi',real=True)
    mix=s.Matrix([[-24*alpha*xi**2,-4/k**2],[-4/k**2,0]])
    ck('mixed Fourier determinant negative',mix.det()+16/k**4)
    ck('zero Fourier h eigenvalue',(s.diff(ell,h,2)+s.diff(ell,hd,2)*xi**2).subs({h:0,V:0})-(6/k**2-24*alpha*xi**2))

    n0,n1,n2,vol=s.symbols('n0 n1 n2 volume',real=True)
    u=s.symbols('u0:4',real=True)
    d=[(n1-n2)**2,(n0-n2)**2,(n0-n1)**2,s.S.Zero]
    Q=s.diag(*d);factor=-4*alpha*vol*h*h
    F=factor*(s.Matrix(u).T*Q*s.Matrix(u))[0]
    lam=s.symbols('lambda',real=True)
    lag=F-factor*lam*(sum(x*x for x in u)-1)
    ck('orientation stationary gradient',sum((s.diff(lag,u[i])-2*factor*(d[i]-lam)*u[i])**2 for i in range(4)))
    ck('orientation constrained Hessian',sum((s.hessian(lag,u)[i,j]-2*factor*(Q[i,j]-lam*(1 if i==j else 0)))**2 for i,j in itertools.product(range(4),repeat=2)))
    rr=s.symbols('r',positive=True)
    ck('Nil orientation spectrum',sum((d[i].subs({n0:0,n1:0,n2:-1/rr})-[1/rr**2,1/rr**2,0,0][i])**2 for i in range(4)))

    return {'Euler_h':str(s.factor(Eh)),'Euler_V':str(s.factor(EV)),
     'mixed_h':str(H),'mixed_V_magnitude':str(K),'resonant_alpha':str(an),
     'second_variation':str(s.factor(second)),'mixed_Fourier_Hessian':str(mix),
     'orientation_matrix':str(Q),'orientation_tangent_Hessian':str(2*factor*(Q-lam*s.eye(4))),
     'identity_residuals':checks}

def coupled_formulas():
    checks={}
    def ck(name,value):checks[name]=str(c.clean(s.simplify(value)))
    r,L,f=s.symbols('r L f',positive=True)
    w,wp,sp,h,V=s.symbols('w w_prime s_prime h V',real=True)
    A=s.diag(r*f,r/f,r,L);A[3,0]=L*w
    Ap=s.diag(r*f*sp,-r/f*sp,0,0);Ap[3,0]=L*wp
    B=A.inv();D=Ap*B
    C1=c.arr(3,lambda a,b,d:c.clean(D[a,d]*B[3,b]-D[a,b]*B[3,d]))
    LC1,_=c.connection(C1,c.torsion(s.S.Zero,s.S.Zero))
    _,K=c.connection(c.arr(3,lambda *i:s.S.Zero),c.torsion(h,V))
    delta1=c.arr(4,lambda a,b,d,e:c.clean(sum(LC1[a,j,d]*K[j,b,e]+K[a,j,d]*LC1[j,b,e]-LC1[a,j,e]*K[j,b,d]-K[a,j,e]*LC1[j,b,d]+K[a,b,j]*C1[j,d,e] for j in range(4))))
    P1=c.project(delta1)
    ck('frequency-linear torsion difference has zero Weyl projection',c.norm(P1['W']))
    ck('frequency-linear torsion difference has zero four-form projection',c.norm(P1['R3']))
    R2linear=c.norm(P1['R2'])
    ck('frequency-linear extra norm has no vector torsion amplitude',s.diff(R2linear,V))

    # At a cosine zero, w=0 and w''=0. The exact X4 relation fixes the spatial second jets.
    v=s.symbols('v',real=True)
    A0=s.diag(r,r,r,L);Ap0=s.zeros(4);Ap0[3,0]=L*v
    App0=s.diag(-L**2*v**2/(2*r),L**2*v**2/(2*r),0,0)
    B0=A0.inv();Bp0=-B0*Ap0*B0
    C10=c.arr(3,lambda a,b,d:s.simplify(sum(Ap0[a,j]*(B0[3,b]*B0[j,d]-B0[3,d]*B0[j,b]) for j in range(4))))
    C1p=c.arr(3,lambda a,b,d:s.simplify(sum(App0[a,j]*(B0[3,b]*B0[j,d]-B0[3,d]*B0[j,b])+Ap0[a,j]*(Bp0[3,b]*B0[j,d]+B0[3,b]*Bp0[j,d]-Bp0[3,d]*B0[j,b]-B0[3,d]*Bp0[j,b]) for j in range(4))))
    G1,_=c.connection(C10,c.torsion(s.S.Zero,s.S.Zero));G1p,_=c.connection(C1p,c.torsion(s.S.Zero,s.S.Zero))
    Rhighest=c.arr(4,lambda a,b,d,e:c.clean(sum(G1[a,j,d]*G1[j,b,e]-G1[a,j,e]*G1[j,b,d]+G1[a,b,j]*C10[j,d,e] for j in range(4))+B0[3,d]*G1p[a,b,e]-B0[3,e]*G1p[a,b,d]))
    highest=c.project(Rhighest)
    q4=c.norm(highest['W'])
    ck('crossing quartic coefficient',q4-v**4/(3*r**4))
    ck('crossing highest extra nine-component sector zero',c.norm(highest['R2']))
    ck('crossing highest extra four-form zero',c.norm(highest['R3']))

    # Direct symbolic identity for arbitrary static shape/twist on this diagonal family.
    C0=c.structure([4*f*f/r,4/(f*f*r),4/r]);C0[3,1,2]=4*L*w*f/r**2;C0[3,2,1]=-C0[3,1,2]
    LC0,_=c.connection(C0,c.torsion(s.S.Zero,s.S.Zero))
    delta0=c.arr(4,lambda a,b,d,e:c.clean(sum(LC0[a,j,d]*K[j,b,e]+K[a,j,d]*LC0[j,b,e]+K[a,j,d]*K[j,b,e]-LC0[a,j,e]*K[j,b,d]-K[a,j,e]*LC0[j,b,d]-K[a,j,e]*K[j,b,d]+K[a,b,j]*C0[j,d,e] for j in range(4))))
    P0=c.project(delta0)
    ck('static torsion difference has zero Weyl projection',c.norm(P0['W']))
    ck('static extra-nine norm',c.norm(P0['R2'])-s.Rational(8,3)*h*h*V*V)
    ck('static four-form norm',c.norm(P0['R3'])-s.Rational(8,3)*h*h*V*V)
    cross2=2*c.inner(P0['R2'],P1['R2'])
    ck('static and frequency-linear extra-nine sectors orthogonal',cross2)
    ck('static scalar torsion shift',s.trace(c.ricci(delta0))+6*h*h+s.Rational(2,3)*V*V)
    ck('frequency-linear scalar torsion shift zero',s.trace(c.ricci(delta1)))
    epsilon_terms=[(p,s.LeviCivita(*p)) for p in itertools.permutations(range(4))]
    ck('static curvature contribution to complete NY',sum(sign*delta0[p]/2 for p,sign in epsilon_terms)+4*h*V)
    ck('frequency-linear complete NY zero',sum(sign*delta1[p]/2 for p,sign in epsilon_terms))
    F=1-L**2*w*w/(2*r*r)
    X4extra=4*h*h*wp*wp*(r*r+L*L*w*w/2)/(r**4*F**2)
    ck('X4 extra-nine derivative coefficient',R2linear.subs({f:s.sqrt(F),sp:-L*L*w*wp/(2*r*r*F)})-X4extra)

    return {'extra_nine_derivative_norm':str(R2linear),
     'extra_nine_static_derivative_cross':str(cross2),'X4_extra_n2_density':str(X4extra),
     'crossing_quartic_W_norm':str(q4),
     'crossing_curvature_nonzero':{','.join(map(str,p)):str(Rhighest[p]) for p in itertools.product(range(4),repeat=4) if Rhighest[p]!=0},
     'identity_residuals':checks}

def nil_quotient_formulas():
    checks={}
    def ck(name,value):checks[name]=str(c.clean(s.simplify(value)))
    a=s.symbols('a',real=True);r,L=s.symbols('r L',positive=True)
    QD=s.diag(s.exp(a),s.exp(-a));QS=s.Matrix([[s.cosh(a),s.sinh(a)],[s.sinh(a),s.cosh(a)]])
    for name,Q in [('T2',QD),('T3',QS)]:
        ck(name+' determinant',Q.det()-1)
        ck(name+' inverse',sum((s.simplify(Q.inv()[i,j])-Q.subs(a,-a)[i,j])**2 for i,j in itertools.product(range(2),repeat=2)))
    Q00,Q01,Q11,x,w,z,y=s.symbols('Q00 Q01 Q11 x w z y',real=True)
    Q=s.Matrix([[Q00,Q01],[Q01,Q11]])
    G=s.Matrix([[r*r*Q00,r*r*Q01,0,0],[r*r*Q01,r*r*(Q11+x*x),-r*r*x,0],[0,-r*r*x,r*r,0],[0,0,0,L*L]])
    Gi=G.inv()
    ck('coordinate determinant',G.det()-L*L*r**6*Q.det())
    ck('horizontal inverse block',sum((Gi[i,j]-Q.inv()[i,j]/r**2)**2 for i,j in itertools.product(range(2),repeat=2)))
    for j in range(2):ck('closed coordinate form coclosed '+str(j),sum(s.diff(Gi[i,j],coord) for i,coord in enumerate([x,w,z,y])))
    R=c.curvature(c.connection(c.structure([0,0,-1/r]),c.torsion(0,0))[0],c.structure([0,0,-1/r]))
    ck('LC Ricci distinguished distributions',sum((c.ricci(R)[i,j]-s.diag(-1/(2*r*r),-1/(2*r*r),1/(2*r*r),0)[i,j])**2 for i,j in itertools.product(range(4),repeat=2)))

    # Integer-lattice automorphisms in normalized (X,Y,Z) coordinates, theta=dZ-X dY.
    X,Y,Z=s.symbols('X Y Z');coordinates=s.Matrix([X,Y,Z]);theta=s.Matrix([0,-X,1])
    rotation=s.Matrix([Y,-X,Z-X*Y]);reflection=s.Matrix([-X,Y,-Z])
    for name,F,sign in [('rotation',rotation,1),('reflection',reflection,-1)]:
        pulled=F.jacobian(coordinates).T*s.Matrix([0,-F[0],1])
        ck(name+' central one-form pullback',sum((pulled-sign*theta)[i]**2 for i in range(3)))
    P=s.Matrix([[0,1],[-1,0]]);J=s.diag(-1,1)
    ck('rotation flips T2 parameter',sum((P.T*QD*P-QD.subs(a,-a))[i,j]**2 for i,j in itertools.product(range(2),repeat=2)))
    ck('reflection flips T3 parameter',sum((J.T*QS*J-QS.subs(a,-a))[i,j]**2 for i,j in itertools.product(range(2),repeat=2)))
    xx,yy,zz=s.symbols('xx yy zz')
    def product(p,q):return s.Matrix([p[0]+q[0],p[1]+q[1],p[2]+q[2]+p[0]*q[1]])
    def apply(F,p):return F.subs({X:p[0],Y:p[1],Z:p[2]},simultaneous=True)
    other=s.Matrix([xx,yy,zz])
    for name,F in [('rotation',rotation),('reflection',reflection)]:
        ck(name+' group automorphism',sum((apply(F,product(coordinates,other))-product(F,apply(F,other)))[i]**2 for i in range(3)))
    ck('T3 threshold',s.cosh(s.log(3)/2)-2*s.exp(-s.log(3)/2))
    return {'Q_T2':str(QD),'Q_T3':str(QS),'horizontal_inverse_metric':str(Q.inv()/r**2),
     'sigma_T2':'exp(-Abs(a))','sigma_T3':'Min(cosh(a),2*exp(-Abs(a)))','a':'sqrt(2)*q',
     'identification_T2':str(rotation),'identification_T3':str(reflection),
     'identity_residuals':checks}

def derive_extensions():
    result={'schema':'paper03ec.extensions.v1','stationary':stationary_formulas(),
            'coupled':coupled_formulas(),'nil_quotient':nil_quotient_formulas()}
    checks={group+'.'+name:value for group in ['stationary','coupled','nil_quotient']
            for name,value in result[group]['identity_residuals'].items()}
    failed={name:value for name,value in checks.items() if value!='0'}
    if failed:raise ArithmeticError(f'Exact identity residuals: {failed}')
    result['identity_residuals']=checks
    return result
