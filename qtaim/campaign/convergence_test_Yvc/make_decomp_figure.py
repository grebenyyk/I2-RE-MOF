import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
DED3={("Gd","st"):-62.90,("Tb","st"):-63.12,("Y","st"):-64.01,
      ("Gd","act"):-62.28,("Tb","act"):-63.51,("Y","act"):-64.48,
      ("Gd","vc"):-84.50,("Tb","vc"):-81.38,("Y","vc"):-83.46}
DAD={("Gd","st"):-151.22,("Tb","st"):-131.74,("Y","st"):-130.08,
     ("Gd","act"):-89.12,("Tb","act"):-77.66,("Y","act"):-69.64,
     ("Gd","vc"):-117.18,("Tb","vc"):-113.71,("Y","vc"):-111.98}
metals=["Y","Tb","Gd"]; variants=[("st","solvent-templated"),("act","activated"),("vc","variable-cell")]
vols={"st":1432,"act":1432,"vc":1242}
fig,axes=plt.subplots(1,3,figsize=(13,6.5),sharey=True)
for ax,(v,vname) in zip(axes,variants):
    x=np.arange(len(metals))
    disp=np.array([DED3[(m,v)] for m in metals])
    elec=np.array([DAD[(m,v)]-DED3[(m,v)] for m in metals])
    ax.bar(x,disp,width=0.55,color='#9ecae1',edgecolor='#3182bd',label='dispersion (D3)  — pore size / wall contact',zorder=3)
    ax.bar(x,elec,width=0.55,bottom=disp,color='#fdae6b',edgecolor='#e6550d',label='electronic (charge transfer)  — I$_2$ oxidation',zorder=3)
    for i,m in enumerate(metals):
        ax.text(i,disp[i]/2,f'{disp[i]:.0f}',ha='center',va='center',fontsize=9,color='#3182bd',fontweight='bold')
        ax.text(i,(disp[i]+elec[i]/2),f'{elec[i]:.0f}',ha='center',va='center',fontsize=9,color='#e6550d',fontweight='bold')
        ax.text(i,(disp[i]+elec[i])-4,f'$\\Sigma$={disp[i]+elec[i]:.0f}',ha='center',va='top',fontsize=9.5,color='k')
    ax.set_xticks(x); ax.set_xticklabels([f'{m}-MOF' for m in metals],fontsize=12)
    ax.set_title(f'{vname}  (V={vols[v]} Å$^3$)',fontsize=12)
    ax.axhline(0,color='k',lw=0.8)
    ax.grid(axis='y',alpha=0.25,zorder=0)
axes[0].set_ylabel(r'$\Delta E_{\rm ads}$ components  (kJ mol$^{-1}$ per I$_2$)',fontsize=12)
axes[0].legend(loc='lower left',fontsize=9.5,framealpha=0.95)
fig.suptitle('Binding-energy decomposition: dispersion (pore size) vs electronic (charge transfer)\n'
             'st & act share an identical cell (1432 Å³) yet differ 2× in binding — the electronic term is the discriminator',
             fontsize=12.5)
plt.tight_layout(rect=[0,0,1,0.93])
out='/workspace/dimitry/calc/qe/i2_mof/HSE_final_tables/landscape_3variants/binding_decomposition.png'
fig.savefig(out,dpi=150,bbox_inches='tight'); fig.savefig('binding_decomposition.png',dpi=150,bbox_inches='tight')
print('saved',out)
