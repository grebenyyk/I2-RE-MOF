import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
V=['st','act','vc']; x=np.arange(3)
anchor =[-132.23,-61.74,-123.40]
assembly=[+2.15,-7.90,-31.36]
deform =[0,0,42.77]
ct=[0.741,0.262,1.392]
elec=[-66.07,-5.16,-28.52]
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,6))
# Panel A: Born-Haber stacked
b1=ax1.bar(x,assembly,width=0.55,color='#fdae6b',edgecolor='#e6550d',label='chain assembly ($\\Delta E_{assembly}$)',zorder=3)
b2=ax1.bar(x,anchor,width=0.55,bottom=assembly,color='#9ecae1',edgecolor='#3182bd',label='anchoring ($\\Delta E_{anchor}$)',zorder=3)
b3=ax1.bar(x,deform,width=0.55,bottom=[assembly[i]+anchor[i] for i in range(3)],color='#d62728',edgecolor='#962b1c',label='framework strain ($\\Delta E_{def}$)',zorder=3)
for i in range(3):
    ax1.text(i,assembly[i]/2,f'{assembly[i]:+.0f}',ha='center',va='center',fontsize=9,color='#7a3b08',fontweight='bold')
    ax1.text(i,assembly[i]+anchor[i]/2,f'{anchor[i]:.0f}',ha='center',va='center',fontsize=9,color='#1f4e79',fontweight='bold')
    if deform[i]>0: ax1.text(i,assembly[i]+anchor[i]+deform[i]/2,f'+{deform[i]:.0f}',ha='center',va='center',fontsize=9,color='white',fontweight='bold')
    tot=assembly[i]+anchor[i]+deform[i]
    ax1.text(i,tot-3,f'$\\Sigma$={tot:.0f}',ha='center',va='top',fontsize=10,fontweight='bold')
ax1.axhline(0,color='k',lw=0.8); ax1.set_xticks(x); ax1.set_xticklabels([f'Y-MOF ({v})' for v in V],fontsize=12)
ax1.set_ylabel('kJ mol$^{-1}$ per I$_2$',fontsize=12); ax1.set_title('Born–Haber: same $\\Delta E_{ads}$, three terms',fontsize=12)
ax1.legend(loc='lower left',fontsize=9.5,framealpha=0.95); ax1.grid(axis='y',alpha=0.25,zorder=0)
# Panel B: the paradox — anchoring/elec vs charge transfer
w=0.32
ax2.bar(x-w/2,anchor,width=w,color='#3182bd',label='$\\Delta E_{anchor}$ (left)',zorder=3)
ax2.bar(x+w/2,elec,width=w,color='#9ecae1',label='$\\Delta E_{elec}$ (left)',zorder=3)
ax2.set_ylabel('$\\Delta E_{anchor},\\ \\Delta E_{elec}$ (kJ mol$^{-1}$ per I$_2$)',color='#3182bd',fontsize=11.5)
ax2.tick_params(axis='y',labelcolor='#3182bd')
ax2.axhline(0,color='#3182bd',lw=0.6)
axr=ax2.twinx()
axr.plot(x,ct,'o-',color='#d62728',lw=2.2,ms=10,label='charge on I (right)',zorder=5)
for i in range(3): axr.annotate(f'{ct[i]:.2f}',(x[i],ct[i]),textcoords='offset points',xytext=(8,6),color='#d62728',fontsize=10,fontweight='bold')
axr.set_ylabel('charge transferred to I (e)',color='#d62728',fontsize=11.5); axr.tick_params(axis='y',labelcolor='#d62728')
ax2.set_xticks(x); ax2.set_xticklabels([f'Y-MOF ({v})' for v in V],fontsize=12)
ax2.set_title('Anchoring/elec $\\uparrow$ st>vc>act   but   charge transfer $\\uparrow$ vc>st>act',fontsize=11.5)
h1,l1=ax2.get_legend_handles_labels(); h2,l2=axr.get_legend_handles_labels()
ax2.legend(h1+h2,l1+l2,loc='upper right',fontsize=9.5,framealpha=0.95)
ax2.grid(axis='y',alpha=0.2,zorder=0)
fig.suptitle('Why vc transfers the most charge but st binds strongest: anchoring (not charge-shuffling) sets the electronic energy',
             fontsize=12.5)
plt.tight_layout(rect=[0,0,1,0.94])
out='/workspace/dimitry/calc/qe/i2_mof/HSE_final_tables/landscape_3variants/anchoring_vs_chargetransfer.png'
fig.savefig(out,dpi=150,bbox_inches='tight'); fig.savefig('anchoring_vs_chargetransfer.png',dpi=150,bbox_inches='tight')
print('saved',out)
