import numpy as np, pandas as pd
df=pd.read_parquet("sim/sim_results_records.parquet")
tau=float(open("sim/sim_tau.txt").read().split("=")[-1])
t=df[df.set=="test"].copy(); t["_ts"]=t.true_scale
t["adeq"]=np.where(t.assigned.str.startswith("abstain"),t.assigned,np.where(t.d_min>tau,"abstain_adequacy",t.assigned))
def decomp(g,col):
    n=len(g); a=g[col]; ts=g["_ts"].iloc[0]
    return pd.Series(dict(n=n,correct=(a==ts).mean(),
        to_linear=((a=="linear_probability")&(ts!="linear_probability")).mean(),
        to_standardised=((a=="standardised")&(ts!="standardised")).mean(),
        to_log_odds=((a=="log_odds")&(ts!="log_odds")).mean(),
        abstain_sufficiency=(a=="abstain_sufficiency").mean(),
        abstain_adequacy=(a=="abstain_adequacy").mean()))
keys=["arm","level","p_target","K","p_supplied_as","pipeline","true_scale"]
P=t.groupby(keys).apply(lambda g:decomp(g,"assigned"),include_groups=False).reset_index().assign(rule="primary (forced nearest)")
A=t.groupby(keys).apply(lambda g:decomp(g,"adeq"),include_groups=False).reset_index().assign(rule="secondary (adequacy abstention)")
U=pd.concat([P,A]); U.to_csv("sim/gate_unconditional_decomposition.csv",index=False)
fit=t.groupby(keys).agg(ge5_gws=("n_pass",lambda s:(s>=5).mean()),median_n_pass=("n_pass","median"),
    d_min_median=("d_min","median"),d_min_p95=("d_min",lambda s:s.quantile(.95)),
    rel_sep_median=("rel_sep","median"),rel_sep_p05=("rel_sep",lambda s:s.quantile(.05))).reset_index()
fit.to_csv("sim/gate_fit_selection.csv",index=False)
ds=df.drop_duplicates(["set","cell","rep"])
nc=ds.groupby(["set","cell"]).conv_logit.agg(lambda s:int((~s).sum())).rename("logistic_nonconverged").reset_index()
nc[nc.logistic_nonconverged>0].to_csv("sim/gate_logistic_nonconvergence.csv",index=False)
print("tau",round(tau,4))
print("\n=== CORE, primary rule, pooled over p, by K x pipeline (unconditional)")
c=P[P.arm=="core"].groupby(["K","pipeline"])[["correct","to_linear","to_standardised","to_log_odds","abstain_sufficiency"]].mean()
print(c.round(3).to_string())
print("\n=== CORE worst cell per K (primary): min correct and where")
w=P[P.arm=="core"].sort_values("correct").groupby("K").head(2)[["K","p_target","pipeline","correct","abstain_sufficiency"]]
print(w.round(3).to_string(index=False))
print("\n=== CORE adequacy rule: false abstention (correct scale abstained) by K x pipeline")
ca=A[A.arm=="core"].groupby(["K","pipeline"]).abstain_adequacy.mean().unstack().round(3); print(ca.to_string())
print("\n=== selection: fraction of releases with >=5 GWS variants, core, by K x p (min over pipelines)")
print(fit[fit.arm=="core"].groupby(["K","p_target"]).ge5_gws.min().unstack().round(3).to_string())
print("\n=== logistic non-convergence: cells with any:",len(nc[nc.logistic_nonconverged>0]),"| total data sets",int(nc.logistic_nonconverged.sum()))
print(nc[nc.logistic_nonconverged>0].sort_values("logistic_nonconverged",ascending=False).head(8).to_string(index=False))
