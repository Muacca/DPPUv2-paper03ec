"""Generate the exact results and fixed datasets for the manuscript figures."""
from pathlib import Path
import hashlib,json
import numpy as np
from .model import circle_components,nil_shape_coefficients

def _csv(path,values,header):
    with path.open('w',encoding='utf-8',newline='\n') as stream:
        np.savetxt(stream,values,delimiter=',',header=header,comments='',fmt='%.17g')

def generate(output):
    from .derive import derive
    out=Path(output)
    if out.exists():raise FileExistsError(f'Choose a new output directory: {out}')
    out.mkdir(parents=True)
    exact=derive()
    (out/'exact_results.json').write_text(json.dumps(exact,indent=2)+'\n',encoding='utf-8',newline='\n')
    y=np.linspace(0,1,401);h=.30+.15*np.sin(2*np.pi*y);V=.40+.20*np.cos(2*np.pi*y)
    dh=.30*np.pi*np.cos(2*np.pi*y);dV=-.40*np.pi*np.sin(2*np.pi*y)
    inv=circle_components([2,2,2],h,V,dh,dV);I2=inv['I2'];I3=inv['I3']
    # S3 r=2: R_LC=6 and W_LC=0. L=1, h=eta/r.
    P=inv['P']
    _csv(out/'circle_profile.csv',np.column_stack([y,h,V,dh,dV,I2,I3,I2+I3,-6*dh,P]),'y,h,V,h_dot,V_dot,I2,I3,TF_minus_W,NY,P')
    q=np.linspace(-.35,.35,281)
    z1=np.exp(-4*q/np.sqrt(6));z4=np.cosh(np.sqrt(2)*q)
    U1=z1*z1+z1**4;U23=np.full_like(q,2.);U45=z4*z4+z4**4
    _csv(out/'nil_paths.csv',np.column_stack([q,U1,U23,U23,U45,U45]),'q,U_T1,U_T2,U_T3,U_T4,U_T5')
    gradient,hessian=nil_shape_coefficients(1,1)
    coeff=np.column_stack([np.arange(1,6),gradient,hessian.diagonal()])
    _csv(out/'nil_coefficients.csv',coeff,'basis,first_derivative,second_derivative')
    from .extensions import derive_extensions
    from .coupled import integrate,action_polynomial
    extended=derive_extensions()
    (out/'extension_results.json').write_text(json.dumps(extended,indent=2)+'\n',encoding='utf-8',newline='\n')
    parameters={'amplitude':.5,'r':3.,'L':1.,'eta':1.,'V':1.}
    mean=integrate(512,**parameters)
    payload={'parameters':parameters,'kappa':1.,'quadrature_samples':512,'mean_invariants':mean,
             'meaning':'ascending frequency coefficients before weighting by -alpha times volume; integer frequencies at fixed amplitude'}
    (out/'coupled_coefficients.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
    n=np.arange(1,17);volume=np.pi**2*27/4;round_lc=-4.5*np.pi**2
    lc=[np.polynomial.polynomial.polyval(n,action_polynomial(mean,alpha=a)) for a in [-1e-4,0.,1e-4]]
    tf=np.polynomial.polynomial.polyval(n,action_polynomial(mean,model='EC_TRACEFREE',alpha=-1e-4))
    constant=np.full(n.shape,1e-4*volume*16/27)
    gradient=1e-4*volume*(mean['EC_TRACEFREE'][2]-mean['LC'][2])*n*n
    _csv(out/'coupled_profiles.csv',np.column_stack([n,*lc,tf,constant,gradient,constant+gradient,np.full(n.shape,round_lc)]),
         'n,S_LC_negative,S_LC_zero,S_LC_positive,S_TF_negative,shift_constant,shift_gradient,shift_total,round_LC')
    files=[{'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(out.iterdir())]
    manifest={'schema':'paper03ec.data.v2','paper_version':'v4','circle':{'topology':'S3','r':2,'L':1,'profile':'h=0.30+0.15 sin(2 pi y); V=0.40+0.20 cos(2 pi y)','direction':'z=L*y','boundary':'smooth periodic'},'nil':{'A':1,'B':1,'chart':'M=exp(q T_A), tr(T_A T_B)=delta_AB','meaning':'fixed-r,L path coefficients; not a stationary background or particle masses'},'coupled':{'parameters':parameters,'kappa':1,'alphas':[-1e-4,0,1e-4],'profile':'w_n=A cos(2 pi n y), f_n^2=1-L^2 w_n^2/(2r^2)','quadrature_samples':512,'frequency_range':[1,16],'meaning':'fixed-amplitude off-shell sequence; no global lower bound or physical stability claim'},'files':files}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    return manifest

def load(source):
    root=Path(source);manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    if manifest['schema'] not in ('paper03ec.data.v1','paper03ec.data.v2'):raise ValueError('Unsupported dataset schema')
    for row in manifest['files']:
        p=(root/row['name']).resolve()
        if not p.is_relative_to(root.resolve()) or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Dataset integrity mismatch')
    names=['circle_profile','nil_paths','nil_coefficients']
    if manifest['schema']=='paper03ec.data.v2':names.append('coupled_profiles')
    return {name:np.genfromtxt(root/(name+'.csv'),delimiter=',',names=True) for name in names}
