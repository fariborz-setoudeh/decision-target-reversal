#!/usr/bin/env python3
"""Decision-Target Reversal: paper figures and Monte Carlo precision.
Run AFTER the simulation and real-data R scripts, from the repository root:
    python 03_reproduce_paper_figures.py .
Install Python dependencies once:
    python -m pip install numpy pandas scipy matplotlib rdata

No new samples or model fits are generated. Reads saved replication estimates,
reconstructs analytic paired risks, reproduces the 2000-resample precision
calculations and second-moment comparisons, and draws Figures 1, 2 and S1.
The bootstrap uses NumPy PCG64(seed=20260930), configurations in the recorded
order (symmetric n=500,1000,2000, then unequal n=2000), batches of 100, and
multinomial counts for paired replication resampling. Do not replace this RNG
with R's RNG if exact reproduction of the published Monte Carlo summaries is
required. Percentile summaries are conditional on recovery of two crossings.

Outputs: paper_figures/*.pdf, *.png, bootstrap_precision.csv,
second_moment_checks.csv, risk_summary_verified.csv.
--figures-only skips bootstrap resampling; all figures are still reconstructed
from saved estimates. This is useful for changing only presentation.
"""
from pathlib import Path
import argparse, hashlib, warnings
import numpy as np
import pandas as pd
from scipy.special import ndtr, expit
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import rdata

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("project_dir", nargs="?", default=".")
parser.add_argument("--figures-only", action="store_true")
args=parser.parse_args()
ROOT=Path(args.project_dir).resolve()
OUT=ROOT/"paper_figures";OUT.mkdir(parents=True,exist_ok=True)
SRC=ROOT/"decision_target_reversal_results"

def locate(directory, filename):
    direct=directory/filename
    if direct.exists(): return direct
    found=sorted(directory.rglob(filename))
    if not found: raise FileNotFoundError(f"Missing {filename} under {directory}. Run the corresponding R script first.")
    signatures={hashlib.sha256(p.read_bytes()).hexdigest() for p in found}
    if len(signatures)!=1: raise ValueError(f"Different copies of {filename}; retain one intended result set.")
    return found[0]

sym=pd.read_csv(locate(SRC/"finite_sample_simulation","all_replicate_parameters_5000.csv"))
uneq=pd.read_csv(locate(SRC/"unequal_prior_robustness","unequal_prior_raw_results_CORRECTED.csv"))
stored=pd.read_csv(locate(SRC/"finite_sample_simulation","finite_sample_risk_summary_5000.csv"))
us=pd.read_csv(locate(SRC/"unequal_prior_robustness","unequal_prior_risk_summary_CORRECTED.csv"))

def diffmatrix(cc,pc,c,p):
    c=np.atleast_1d(c);lc=np.log(c)
    tc=(lc[None,:]-cc[:,0,None])/cc[:,1,None]
    tp=(lc[None,:]-pc[:,0,None])/pc[:,1,None]
    return p*(ndtr(tp-.75)-ndtr(tc-.75))+c[None,:]*(1-p)*(ndtr(-tp-.75)-ndtr(-tc-.75))

def crossings(x,y):
    ids=np.flatnonzero(y[:-1]*y[1:]<0)
    roots=list(x[ids]-y[ids]*(x[ids+1]-x[ids])/(y[ids+1]-y[ids]))
    # Include exact isolated zeros if the grid happens to contain one.
    for j in range(1,len(y)-1):
        if y[j]==0 and y[j-1]*y[j+1]<0:roots.append(x[j])
    return np.exp(np.sort(np.unique(roots)))

rng=np.random.default_rng(20260930);boot=2000
checks=[];precision=[];riskrows=[]
configs=[(f"symmetric_{n}",sym[sym.n==n].sort_values("rep"),.5,n,
          ["cc_beta0","cc_Delta","pc_beta0","pc_Delta"]) for n in [500,1000,2000]]+[
    ("unequal_2000",uneq.sort_values("rep"),.3,2000,["beta0_CC","Delta_CC","beta0_PC","Delta_PC"])]
for name,df,p,n,cols in configs:
    B=len(df);expected=5000 if p==.5 else 3000
    if B!=expected or df.rep.nunique()!=B or set(df.rep)!=set(range(1,B+1)):
        raise ValueError(f"Incomplete or duplicate replication IDs: {name}")
    success=df.paired_success if p==.5 else df.success
    if not success.all():raise ValueError(f"{name} contains unsuccessful fits; do not silently condition on success")
    cc=df[cols[:2]].to_numpy();pc=df[cols[2:]].to_numpy();true=np.array([np.log(p/(1-p)),1.5])
    if not np.isfinite(np.r_[cc.ravel(),pc.ravel()]).all() or np.any(cc[:,1]<=0) or np.any(pc[:,1]<=0):
        raise ValueError(f"Invalid parameter estimates: {name}")
    ec=cc-true;ep=pc-true
    D=n*(np.cov(pc.T)-np.cov(cc.T))
    M=n*(ep.T@ep-ec.T@ec)/B
    pop=np.array([[-.471312363222559,0],[0,.447127973747152]]) if p==.5 else np.array([
        [-.095231683321602,-.328968287243148],[-.328968287243148,.199038600804896]])
    orig=(stored[stored.n==n] if p==.5 else us).sort_values("c")
    costs=orig.c.to_numpy();v=diffmatrix(cc,pc,costs,p)
    means=v.mean(axis=0);ses=v.std(axis=0,ddof=1)/np.sqrt(B)
    error=np.max(abs(means-orig.mean_diff.to_numpy()))
    if error>1e-12:raise ValueError(f"Saved mean risks do not match parameter estimates: {name}; discrepancy={error}")
    t=(np.log(costs)-true[0])/1.5
    pref=p*np.exp(-.5*(t-.75)**2)/np.sqrt(2*np.pi)/3/n
    approx=lambda A:pref*(A[0,0]+2*t*A[0,1]+t*t*A[1,1])
    approxD,approxM,approxP=approx(D),approx(M),approx(pop)
    checks.append(dict(configuration=name,n=n,replications=B,
        relative_second_moment_change_percent=100*np.linalg.norm(M-D)/np.linalg.norm(D),
        RMSE_centered_covariance=np.sqrt(np.mean((means-approxD)**2)),
        RMSE_second_moment=np.sqrt(np.mean((means-approxM)**2)),
        RMSE_population=np.sqrt(np.mean((means-approxP)**2)),
        maximum_risk_reproduction_error=error))
    for j,c in enumerate(costs):
        riskrows.append(dict(configuration=name,n=n,B=B,c=c,log_c=np.log(c),mean=means[j],
            MCSE=ses[j],lower=means[j]-1.96*ses[j],upper=means[j]+1.96*ses[j],
            population_theory=approxP[j],covariance_theory=approxD[j]))
    if not args.figures_only:
        samples=[]
        for first in range(0,boot,100):
            count=min(100,boot-first)
            weights=rng.multinomial(B,np.full(B,1/B),size=count).astype(float)/B
            for y in weights@v:
                roots=crossings(np.log(costs),y)
                samples.append(roots if len(roots)==2 else [np.nan,np.nan])
        samples=np.asarray(samples);recovered=np.isfinite(samples).all(axis=1)
        ranges=np.quantile(samples[recovered],[.025,.975],axis=0) if recovered.any() else np.full((2,2),np.nan)
        precision.append(dict(configuration=name,resamples=boot,two_crossings_fraction=recovered.mean(),
            lower_crossing_2p5=ranges[0,0],lower_crossing_97p5=ranges[1,0],
            upper_crossing_2p5=ranges[0,1],upper_crossing_97p5=ranges[1,1]))
        print(name,precision[-1],flush=True)

pd.DataFrame(checks).to_csv(OUT/"second_moment_checks.csv",index=False)
if precision:pd.DataFrame(precision).to_csv(OUT/"bootstrap_precision.csv",index=False)
risk=pd.DataFrame(riskrows);risk.to_csv(OUT/"risk_summary_verified.csv",index=False)

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"pdf.fonttype":42,
    "axes.spines.top":False,"axes.spines.right":False})
def simulation_figure(configuration, low, high, basename):
    d=risk[risk.configuration==configuration].sort_values("log_c");n=int(d.n.iloc[0]);B=int(d.B.iloc[0])
    x=d.log_c.to_numpy()
    fig,ax=plt.subplots(figsize=(7.4,4.8))
    ax.axvspan(np.log(low),np.log(high),color="grey",alpha=.065,lw=0)
    for boundary in [low,high]:ax.axvline(np.log(boundary),color="grey",lw=.8,ls=":")
    ax.fill_between(x,n*d.lower.to_numpy(),n*d.upper.to_numpy(),color="#56B4E9",alpha=.23,
        label="95% pointwise Monte Carlo interval")
    ax.plot(x,n*d["mean"].to_numpy(),color="#0072B2",lw=1.9,label="Monte Carlo mean")
    ax.plot(x,n*d.covariance_theory.to_numpy(),color="#D55E00",lw=1.5,ls="--",
        label="Quadratic comparison: empirical covariance")
    ax.plot(x,n*d.population_theory.to_numpy(),color="#009E73",lw=1.6,ls=":",
        label="Population first-order theory")
    ax.axhline(0,color="grey",lw=.8,ls="--")
    ax.set_xlabel(r"$\log(C_{12}/C_{21})$");ax.set_ylabel(r"$n\,\widehat{\Delta R}_n(c)$")
    ax.set_title(f"n = {n:,}; {B:,} paired replications",fontsize=12,pad=13)
    ax.legend(loc="upper left",frameon=False,fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT/(basename+".pdf"));fig.savefig(OUT/(basename+".png"),dpi=600,facecolor="white")
    plt.close(fig)
simulation_figure("symmetric_2000",.214374234287,4.66473969377,"Figure_Main_With_MC_Uncertainty")
simulation_figure("unequal_2000",.347964,75.145039,"Figure_Unequal_With_MC_Uncertainty")

# Figure 2: use only the checked common-covariance reanalysis.
gastro=ROOT/"review_gastro_common_covariance"
pred=pd.read_csv(gastro/"predictions.csv").sort_values("fold")
status=pd.read_csv(gastro/"fold_status.csv").sort_values("fold")
if len(pred)!=76 or list(pred.fold)!=list(range(1,77)) or len(status)!=76 or not status.paired_pass.all():
    raise ValueError("All 76 paired real-data fits must pass before plotting")
if not pred.paired_pass.all():raise ValueError("Unvalidated real-data predictions")
with warnings.catch_warnings():
    # SessionInfo classes are irrelevant to extracting the saved input matrix.
    warnings.simplefilter("ignore",UserWarning)
    saved=rdata.read_rds(gastro/"common_covariance_results.rds")
X=np.asarray(saved["input"]["X"])
truth=pred.truth.to_numpy();observed=pred.label_available.to_numpy(dtype=bool)
if X.shape!=(76,4) or not np.array_equal(truth,np.asarray(saved["input"]["truth"])):
    raise ValueError("Real-data prediction order does not match the saved input")
if not np.array_equal(observed,~np.isnan(np.asarray(saved["input"]["zm"],dtype=float))):
    raise ValueError("Label mask differs from the saved input")
lc=pred.CC_logodds.to_numpy();lp=pred.FULL_validated_logodds.to_numpy()
if not np.isfinite(np.r_[lc,lp]).all():raise ValueError("Nonfinite held-out log-odds")
def coefficients_log(x):
    cc=lc>x;full=lp>x
    return (np.sum(~full&(truth==1))-np.sum(~cc&(truth==1)),
            np.sum(full&(truth==2))-np.sum(cc&(truth==2)))
def risks_log(x):
    cc=lc>x;full=lp>x;c=np.exp(x)
    return ((np.sum(~cc&(truth==1))+c*np.sum(cc&(truth==2)))/76,
            (np.sum(~full&(truth==1))+c*np.sum(full&(truth==2)))/76)
edges=np.unique(np.r_[-3,np.r_[lc,lp][abs(np.r_[lc,lp])<3],3]);pieces=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    a,b=coefficients_log((lo+hi)/2);root=-a/b if b else np.nan
    logroot=np.log(root) if np.isfinite(root) and root>0 else np.nan
    cuts=[lo]+([logroot] if lo<logroot<hi else [])+[hi]
    for x,y in zip(cuts[:-1],cuts[1:]):
        pieces.append(dict(left=np.exp(x),right=np.exp(y),log_right=y,A=a,B=b,
                          sign=int(np.sign(a+b*np.exp((x+y)/2)))))
trans=[]
for x,y in zip(pieces[:-1],pieces[1:]):
    if x["sign"]!=y["sign"]:trans.append(dict(c=x["right"]))
print("Equal-cost real-data risks:",risks_log(0),flush=True)
# Log entropy from held-out FULL predictions, stable at extreme log odds.
t=np.abs(lp);H=np.logaddexp(0,-t)+t*expit(-t)
for group in [True,False]:print('entropy group',group,'n',np.sum(observed==group),'mean',H[observed==group].mean(),'median',np.median(H[observed==group]))
Xs=(X-X.mean(0))/X.std(0,ddof=1);u,sv,vt=np.linalg.svd(Xs,full_matrices=False);pc=u*sv;expl=100*sv**2/sum(sv**2)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'axes.linewidth':.8,'axes.edgecolor':'#363B40','text.color':'#252A30','axes.labelcolor':'#252A30','pdf.fonttype':42})
fig,(a,b,c)=plt.subplots(1,3,figsize=(11.6,4.5),gridspec_kw={'width_ratios':[1.12,.85,1.62]})
fig.subplots_adjust(left=.061,right=.988,bottom=.215,top=.935,wspace=.37)
colors={1:'#18AEB5',2:'#E76F61'}
for cl in [2,1]:
 for available,marker,size in [(True,'o',36),(False,'^',48)]:
  keep=(truth==cl)&(observed==available);a.scatter(pc[keep,0],pc[keep,1],color=colors[cl],marker=marker,s=size,edgecolors='white',linewidths=.45,alpha=.93,zorder=3)
a.set_xlabel(f'PC1 ({expl[0]:.1f}%)');a.set_ylabel(f'PC2 ({expl[1]:.1f}%)');a.margins(.08)
rng=np.random.default_rng(20260930)
for position,keep,color in zip([0,1],[observed,~observed],['#526D82','#D99A3D']):
 h=H[keep];grid=np.linspace(0,np.log(2),300);den=gaussian_kde(h)(grid);width=.33*den/den.max()
 b.fill_betweenx(grid,position-width,position+width,color=color,alpha=.18,linewidth=0)
 b.plot(position-width,grid,color=color,lw=1);b.plot(position+width,grid,color=color,lw=1)
 b.scatter(position+rng.uniform(-.055,.055,len(h)),h,s=12,color=color,alpha=.72,lw=0,zorder=3)
 b.boxplot([h],positions=[position],widths=.16,showfliers=False,patch_artist=True,boxprops={'facecolor':'white','edgecolor':color,'alpha':.72},medianprops={'color':color,'linewidth':1.4},whiskerprops={'color':color},capprops={'color':color})
b.set_xticks([0,1],['Available','Masked']);b.set_xlabel('Label availability');b.set_ylabel('Held-out FULL posterior entropy');b.set_xlim(-.48,1.48);b.set_ylim(-.025,.725);b.set_yticks([0,.2,.4,.6])
def plot_pieces(ax,limits,shade=True):
 ax.axhline(0,color='#4E545A',linestyle=(0,(4,3)),lw=.75)
 previous=None
 for s in pieces:
  left,right=s['left'],s['right'];grid=np.geomspace(left,right,max(4,int(80*np.log(right/left))+2));yy=(s['A']+s['B']*grid)/76
  if shade:ax.axvspan(np.log(left),np.log(right),color={-1:'#2A8C82',0:'#92999F',1:'#C76D55'}[s['sign']],alpha=.075 if s['sign'] else .12,lw=0)
  ax.plot(np.log(grid),yy,color='#173F4F',lw=1.45,zorder=4)
  if previous is not None and previous!=yy[0]:ax.plot([np.log(left),np.log(left)],[previous,yy[0]],color='#173F4F',lw=.75,alpha=.8)
  previous=yy[-1]
 ax.set_xlim(limits)
plot_pieces(c,[-3.03,3.03]);c.set_ylim(-.13,.78);c.set_xticks(range(-3,4));c.set_yticks([-.1,0,.2,.4,.6]);c.set_xlabel(r'$\log(c),\quad c=C_{12}/C_{21}$');c.set_ylabel(r'$\widehat R_{\mathrm{FULL}}(c)-\widehat R_{\mathrm{CC}}(c)$')
c.scatter([0],[-3/76],s=27,facecolor='white',edgecolor='#173F4F',lw=1.1,zorder=5)
c.annotate(r'$c=1$',xy=(0,-3/76),xytext=(.13,.035),fontsize=8,arrowprops={'arrowstyle':'-','lw':.65,'color':'#627079'},color='#46545D')
xhi=np.log(trans[-1]['c']);c.axvline(xhi,color='#747B80',linestyle=(0,(1.5,2.5)),lw=.85)
c.text(xhi-.05,.61,r'$c=2.463$',rotation=90,ha='right',va='center',fontsize=8,color='#62696E')
c.text(2.03,.70,'CC favourable',ha='center',fontsize=8,color='#B75D49');c.text(-.19,-.112,'FULL favourable',ha='center',fontsize=8,color='#237B72')
ins=c.inset_axes([.08,.52,.44,.39]);plot_pieces(ins,[-2.14,-1.18]);ins.set_ylim(-.014,.014);ins.set_xticks([-2,-1.6,-1.2]);ins.set_yticks([-.01,0,.01]);ins.tick_params(labelsize=6.5,length=2);ins.set_title('Low-cost detail',fontsize=7.5,pad=3);ins.set_facecolor('white')
for ax,tag in zip([a,b,c],['a','b','c']):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(direction='out',length=3,width=.7);ax.text(-.16,1.015,tag,transform=ax.transAxes,fontweight='bold',fontsize=12,ha='left',va='bottom')
legend=[Line2D([],[],ls='none',label='Histological truth'),Line2D([],[],marker='o',ls='none',color=colors[2],markersize=6,label='No resection'),Line2D([],[],marker='o',ls='none',color=colors[1],markersize=6,label='Resection'),Line2D([],[],ls='none',label='Label availability'),Line2D([],[],marker='o',ls='none',color='#454B50',markersize=5.5,label='Available'),Line2D([],[],marker='^',ls='none',color='#454B50',markersize=6,label='Masked')]
fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.025),ncol=6,frameon=False,fontsize=9.2,columnspacing=1.15,handlelength=1,handletextpad=.35)
base=OUT/'Figure_Gastro_Common_Covariance_FINAL'
fig.savefig(base.with_suffix('.pdf'),metadata={'Title':'Gastrointestinal application: shared covariance and training-only LOOCV','Subject':'All 76 folds pass numerical checks; held-out entropy and exact empirical risk ordering'})
fig.savefig(base.with_suffix('.png'),dpi=600,facecolor='white')
plt.close(fig)
print('All paper figures saved in',OUT)
