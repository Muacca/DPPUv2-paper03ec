"""Render saved data without evaluating the symbolic model."""
from pathlib import Path
import matplotlib.pyplot as plt
from .data import load

NAMES=('curvature_components','nil_shape_paths','coupled_profiles')

def figure(name,source):
    data=load(source)
    with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120}):
        if name==NAMES[0]:
            d=data['circle_profile'];fig,ax=plt.subplots(1,2,figsize=(9.4,3.25),layout='constrained')
            ax[0].plot(d['y'],d['h'],label=r'$h=\eta/r$',color='#195a91')
            ax[0].plot(d['y'],d['V'],label=r'$V$',color='#bb5a1d')
            ax[0].set(xlabel=r'Circle coordinate $y$',ylabel='Profile amplitude',title='(a) Smooth periodic torsion')
            ax[1].plot(d['y'],d['I2'],label=r'$I_2$',color='#195a91')
            ax[1].plot(d['y'],d['I3'],label=r'$I_3$',color='#bb5a1d')
            ax[1].plot(d['y'],d['TF_minus_W'],label=r'$I_{\rm tf}-I_{\rm W}$',color='#363e45',linestyle='--')
            ax[1].set(xlabel=r'Circle coordinate $y$',ylabel='Squared-curvature density',title='(b) Additional curvature components')
            for a in ax:a.legend(frameon=False);a.grid(alpha=.16)
        elif name==NAMES[1]:
            d=data['nil_paths'];c=data['nil_coefficients'];fig,ax=plt.subplots(1,2,figsize=(9.4,3.25),layout='constrained')
            ax[0].plot(d['q'],d['U_T1'],label=r'$T_1$',color='#195a91')
            ax[0].plot(d['q'],d['U_T2'],label=r'$T_2,T_3$',color='#626a73',linestyle='--')
            ax[0].plot(d['q'],d['U_T4'],label=r'$T_4,T_5$',color='#bb5a1d')
            ax[0].set(xlabel=r'Exponential path parameter $q$',ylabel=r'$U/A$',title='(a) Nil shape paths at $A=B$')
            ax[0].legend(frameon=False);ax[0].grid(alpha=.16)
            ax[1].bar(c['basis']-.17,c['first_derivative'],width=.34,label='First derivative',color='#195a91')
            ax[1].bar(c['basis']+.17,c['second_derivative'],width=.34,label='Second derivative',color='#bb5a1d')
            ax[1].axhline(0,color='#777777',linewidth=.7);ax[1].set_xticks(c['basis'],[r'$T_'+str(i)+'$' for i in range(1,6)])
            ax[1].set(ylabel='Coefficient divided by $A$',title='(b) Specified path coefficients')
            ax[1].legend(frameon=False)
        elif name==NAMES[2]:
            if 'coupled_profiles' not in data:raise ValueError('This figure requires the current dataset with coupled profiles')
            d=data['coupled_profiles'];d=d[d['n']<=8]
            fig,ax=plt.subplots(1,2,figsize=(9.4,3.25),layout='constrained')
            for column,label,color in [('S_LC_negative',r'$\alpha=-10^{-4}$','#195a91'),('S_LC_zero',r'$\alpha=0$','#626a73'),('S_LC_positive',r'$\alpha=10^{-4}$','#bb5a1d')]:
                ax[0].plot(d['n'],d[column]-d['round_LC'],'.-',label=label,color=color)
            ax[0].axhline(0,color='#888888',linewidth=.7)
            ax[0].set(xlabel=r'Integer frequency $n$',ylabel=r'$S_{\rm LC}(n)-S_{\rm round}$',title='(a) LC-Weyl / EC-Weyl action')
            ax[1].plot(d['n'],d['shift_total'],'.-',label='Total shift',color='#363e45')
            ax[1].plot(d['n'],d['shift_gradient'],'--',label='Nine-component derivative',color='#195a91')
            ax[1].plot(d['n'],d['shift_constant'],':',label='Constant mixed term',color='#bb5a1d',linewidth=2)
            ax[1].set(xlabel=r'Integer frequency $n$',ylabel=r'$S_{\rm tf}(n)-S_{\rm W}(n)$',title=r'(b) Additional sectors, $\alpha=-10^{-4}$')
            for a in ax:a.legend(frameon=False,fontsize=9);a.grid(alpha=.16);a.set_xticks([1,2,4,6,8])
        else:raise ValueError(f'Unknown figure: {name}')
    return fig

def render(source,output):
    out=Path(output)
    if out.exists():raise FileExistsError(f'Choose a new output directory: {out}')
    out.mkdir(parents=True)
    names=NAMES if 'coupled_profiles' in load(source) else NAMES[:2]
    for name in names:
        fig=figure(name,source);fig.savefig(out/(name+'.png'),dpi=200);plt.close(fig)
    return [out/(name+'.png') for name in names]
