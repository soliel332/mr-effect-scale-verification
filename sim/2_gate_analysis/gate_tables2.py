import numpy as np, pandas as pd
U=pd.read_csv("sim/gate_unconditional_decomposition.csv")
U["wrong"]=U[["to_linear","to_standardised","to_log_odds"]].sum(1)
df=pd.read_parquet("sim/sim_results_records.parquet"); t=df[df.set=="test"]
core=t[t.arm=="core"]; nw=int((~core.assigned.str.startswith("abstain") & (core.assigned!=core.true_scale)).sum())
print("CORE exact count of wrong labels (primary):",nw,"of",len(core),"releases")
arms=U[U.arm.isin(["A4","A5","A6","A7","A8"])]
pv=arms.pivot_table(index=["arm","level","p_supplied_as","p_target","pipeline"],columns="rule",values=["wrong","correct","abstain_adequacy"])
pv.columns=[f"{a}|{'prim' if b.startswith('primary') else 'adeq'}" for a,b in pv.columns]
pv=pv[["correct|prim","wrong|prim","correct|adeq","wrong|adeq","abstain_adequacy|adeq"]].reset_index()
pv.to_csv("sim/gate_adversarial_primary_vs_adequacy.csv",index=False)
sig=pv[(pv["wrong|prim"]>0.05)]
print("\n=== adversarial cells with primary misattribution >0.05 (n=%d): primary wrong vs adequacy wrong / abstain"%len(sig))
print(sig.round(3).to_string(index=False))
print("\n=== adversarial arms: max primary wrong per arm/level (all p, pipelines)")
print(pv.groupby(["arm","level","p_supplied_as"])[["wrong|prim","wrong|adeq","abstain_adequacy|adeq"]].max().round(3).to_string())
# audit coverage
a=pd.read_csv("deliverables/reproduction/effect-scale-mr-qc/results/repository_scale_reanalysis.csv")
b=pd.cut(a.nSNP,[4,9,29,30,10**6],labels=["5-9","10-29","30","31+"])
print("\n=== audit nSNP (41 classifiable):",sorted(a.nSNP.tolist()))
print(b.value_counts().reindex(["5-9","10-29","30","31+"]).to_string())
print("within simulated K range (5-30):",int((a.nSNP<=30).sum()),"| above 30:",int((a.nSNP>30).sum()))
L=0.5*np.log(a.prev*(1-a.prev)); d=np.abs(np.log(a.sigma_obs)[:,None]-np.c_[L,0*L,-L]).min(1)
a["d_min"]=d; tau=float(open("sim/sim_tau.txt").read().split("=")[-1])
ex=a[a.d_min>tau][["id","pref","nSNP","prev","sigma_obs","inferred_recomputed","d_min","metadata_discordant"]]
print("\n=== POST HOC: audit records whose distance exceeds tau=%.4f (would abstain under the adequacy rule): %d"%(tau,len(ex)))
print(ex.sort_values("d_min",ascending=False).round(4).to_string(index=False))
print("discordant records exceeding tau:",int(ex.metadata_discordant.astype(str).isin(["True","1"]).sum()))
a[["id","pref","nSNP","prev","d_min","inferred_recomputed","metadata_discordant"]].to_csv("sim/gate_audit_coverage_posthoc.csv",index=False)
