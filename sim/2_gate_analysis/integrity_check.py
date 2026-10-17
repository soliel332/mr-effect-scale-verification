import os,sys,hashlib,numpy as np,pandas as pd
sys.path.insert(0,"sim"); import sim_scale as ss, run_sim as rs
I={}
T=rs.tasks("all")
exp={f"{s}__{c['cell']}__{a:04d}.parquet":(s,c,a,b) for s,c,a,b,_ in T}
got=set(os.listdir("sim/out"))
I["expected chunks"]=len(exp); I["chunk files present"]=len([g for g in got if g.endswith(".parquet")])
I["missing chunks"]=len(set(exp)-got); I["unexpected files"]=len(got-set(exp)); I["leftover .tmp"]=len([g for g in got if g.endswith(".tmp")])
parts=[]; bad=0
for f,(s,c,a,b) in sorted(exp.items()):
    x=pd.read_parquet(f"sim/out/{f}"); x["_file"]=f; parts.append(x)
    bad+= len(x)!=(b-a)*7*(2 if c.get("design")=="casecontrol" else 1)
df=pd.concat(parts,ignore_index=True)
I["chunks with wrong row count"]=int(bad); I["total release records"]=len(df)
I["duplicate (set,cell,rep,pipeline,p_supplied_as)"]=int(df.duplicated(["set","cell","rep","pipeline","p_supplied_as"]).sum())
keys=[]
for s,c,reps,seq in rs.seed_plan():
    for k in seq.spawn(reps): keys.append((k.entropy,)+tuple(k.spawn_key))
I["replicate seeds in plan"]=len(keys); I["unique replicate seeds"]=len(set(keys))
ds=df.drop_duplicates(["set","cell","rep"])
I["data sets tuning"]=int((ds.set=="tuning").sum()); I["data sets test"]=int((ds.set=="test").sum())
I["logistic non-converged data sets"]=int((~ds.conv_logit).sum()); I["Firth non-converged data sets"]=int((~ds.conv_firth).sum())
I["data sets with any SPA fallback"]=int((ds.spa_ok_frac<1).sum()); I["mean SPA success fraction"]=round(float(ds.spa_ok_frac.mean()),6)
I["non-finite sigma_hat outside sufficiency abstentions"]=int((~np.isfinite(df.sigma_hat)&(df.assigned!="abstain_sufficiency")).sum())
I["sufficiency abstentions (records)"]=int((df.assigned=="abstain_sufficiency").sum())
fz=open("sim/FREEZE.txt").read()
for f in ["sim/SIMULATION_PROTOCOL.md","sim/sim_scale.py","sim/run_sim.py","sim/analyse_sim.py"]:
    I[f"hash matches FREEZE: {os.path.basename(f)}"]=hashlib.sha256(open(f,"rb").read()).hexdigest() in fz
import scipy,pyarrow
I["versions"]=f"numpy {np.__version__} scipy {scipy.__version__} pandas {pd.__version__} pyarrow {pyarrow.__version__}"
# deterministic regeneration: one replicate from every cell of the plan, compared field-for-field
plan=rs.seed_plan(); rng=np.random.default_rng(99); mism=0; ncmp=0; worst=0.0
for s,c,reps,seq in plan:
    kids=seq.spawn(reps); r=int(rng.integers(reps))
    new=pd.DataFrame(rs.one_replicate(s,c,r,kids[r]))
    old=df[(df.set==s)&(df.cell==c["cell"])&(df.rep==r)].drop(columns="_file")
    m=new.merge(old,on=["pipeline","p_supplied_as"],suffixes=("_n","_o"))
    ncmp+=len(m)
    mism+= int((m.assigned_n!=m.assigned_o).sum()) + int((m.n_pass_n!=m.n_pass_o).sum())
    for col in ["sigma_hat","d_min","separation","rel_sep"]:
        a,b=m[f"{col}_n"].to_numpy(float),m[f"{col}_o"].to_numpy(float)
        both=np.isfinite(a)&np.isfinite(b); mism+=int((np.isfinite(a)!=np.isfinite(b)).sum())
        if both.any(): worst=max(worst,float(np.max(np.abs(a[both]-b[both]))))
I["regeneration: cells sampled"]=len(plan); I["regeneration: records compared"]=ncmp
I["regeneration: label / n_pass / finiteness mismatches"]=mism; I["regeneration: max abs diff in continuous fields"]=worst
s=pd.Series(I); s.to_csv("sim/integrity_check_results.csv",header=["value"]); print(s.to_string())
