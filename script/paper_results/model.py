"""Restricted circle formulas; dots mean d/dz, z=L*y.

The metric and spatial shape are fixed and twist/mixing vanish.
The three selectors change only the quadratic invariant; all retain R_EC.
"""
import numpy as np

MODELS={'LC':'LC-Weyl comparison','EC_IRR':'Irreducible EC-Weyl','EC_TRACEFREE':'Ricci-tracefree EC curvature'}

def circle_components(n,h,V,h_dot=0.,V_dot=0.):
    n=np.asarray(n,dtype=float)
    if n.shape!=(3,):raise ValueError('n must contain three constant spatial structure eigenvalues')
    h,V,hd,Vd=np.broadcast_arrays(h,V,h_dot,V_dot)
    lam=np.array([(n[i]**2-(n[(i+1)%3]-n[(i+2)%3])**2)/2 for i in range(3)])
    R0=lam.sum();I0=2*(lam@lam)-2*R0*R0/3
    I2=6*(hd+2*h*V/3)**2;I3=6*(hd-2*h*V/3)**2
    return {'R_LC':R0,'I_LC':I0,'I_W':I0,'I2':I2,'I3':I3,'I_tf':I0+I2+I3,
            'R_EC':R0-6*h*h-2*V*V/3+2*Vd,'NY':-6*hd,
            'P':(2*R0-12*h*h-4*V*V/3)*hd-8*V*Vd*h/3}

def circle_density(model,n,h,V,h_dot=0.,V_dot=0.,*,kappa=1.,alpha=0.,theta=0.,bulk=False):
    """Return S/vol density, excluding c_M L r^3.

    bulk=True removes -V_dot/kappa^2 + 6 theta h_dot, whose integral is
    zero for smooth periodic profiles and constant geometry.
    """
    if model not in MODELS:raise ValueError('model must be one of '+', '.join(MODELS))
    if kappa<=0:raise ValueError('kappa must be positive')
    inv=circle_components(n,h,V,h_dot,V_dot)
    key={'LC':'I_LC','EC_IRR':'I_W','EC_TRACEFREE':'I_tf'}[model]
    value=-inv['R_EC']/(2*kappa*kappa)-alpha*inv[key]-theta*inv['NY']
    return value+np.asarray(V_dot)/(kappa*kappa)-6*theta*np.asarray(h_dot) if bulk else value

def nil_shape_coefficients(A,B):
    """Gradient and Hessian at q=0 for the normalized exponential shape chart."""
    return np.array([-4*(2*A+4*B)/np.sqrt(6),0,0,0,0]),np.diag([(32*A+128*B)/3,0,0,4*A+8*B,4*A+8*B])
