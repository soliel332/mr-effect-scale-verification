import numpy as np, pandas as pd
df=pd.read_parquet("sim/sim_results_records.parquet"); t=df[(df.set=="test")]
for lvl,p,pl in [("info0.5",0.5,"P1"),("info0.5",0.5,"P2"),("Nx5.0",0.5,"P1"),("Nx5.0",0.5,"P2"),("Nx0.5",0.5,"P2")]:
    g=t[(t.level==lvl)&(t.p_target==p)&(t.pipeline==pl)]
    ab=g.assigned.str.startswith("abstain"); cls=g[~ab]
    print(f"{lvl} p={p} {pl}: n={len(g)} sufficiency-abstained={ab.mean():.3f} | among classified: wrong={(cls.assigned!=cls.true_scale).mean():.3f} (n={len(cls)}) | unconditional wrong={((~ab)&(g.assigned!=g.true_scale)).mean():.3f} | assigned as: {cls.assigned.value_counts().to_dict()}")
a=pd.read_csv("deliverables/reproduction/effect-scale-mr-qc/results/repository_scale_reanalysis.csv")
L=(0.5*np.log(a.prev*(1-a.prev))).to_numpy(); ls=np.log(a.sigma_obs.to_numpy())
a["d_min"]=np.abs(ls[:,None]-np.c_[L,0*L,-L]).min(1)
tau=float(open("sim/sim_tau.txt").read().split("=")[-1])
ex=a[a.d_min>tau][["id","pref","nSNP","prev","sigma_obs","inferred_recomputed","d_min","metadata_discordant"]]
print(f"\nPOST HOC: audit records with d_min > tau={tau:.4f}: {len(ex)} of 41")
print(ex.sort_values("d_min",ascending=False).round(4).to_string(index=False))
print("discordant among them:",int(ex.metadata_discordant.astype(str).isin(["True","1"]).sum()))
dd=a[a.metadata_discordant.astype(str).isin(["True","1"])]
print("\n11 discordant: d_min range",round(dd.d_min.min(),4),"-",round(dd.d_min.max(),4),"| nSNP",sorted(dd.nSNP.tolist()),"| case fractions",sorted(dd.prev.round(3).tolist()))
# N-factor needed to move each discordant record from linear to standardised, i.e. sigma_std/sigma_obs squared
dd=dd.assign(N_factor_to_look_standardised=(1/dd.sigma_obs)**2)
print("factor by which reported N would have to UNDERSTATE the true N for a truly standardised release to appear as observed:",
      dd.N_factor_to_look_standardised.round(2).sort_values().tolist())
a[["id","pref","nSNP","prev","sigma_obs","d_min","inferred_recomputed","metadata_discordant"]].to_csv("sim/gate_audit_coverage_posthoc.csv",index=False)
