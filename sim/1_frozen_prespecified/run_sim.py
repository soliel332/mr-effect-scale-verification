"""
Runner for SIMULATION_PROTOCOL.md. Resume-safe: each (cell, replicate-chunk) is written atomically
to sim/out/ and skipped if already present. Seeds derive from one SeedSequence, so results do not
depend on chunking or worker count.

    python sim/run_sim.py [--workers 10] [--only tuning|test|all]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import sys, argparse, time
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim_scale as ss

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P_CORE = [0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75]
P_ARM = [0.02, 0.10, 0.35, 0.50]
N = 20000


def build_cells():
    core = [dict(cell=f"core_p{p}_K{K}", arm="core", level="-", N=N, K=K, p=p)
            for K in (5, 10, 30) for p in P_CORE]
    arms = []
    def add(arm, level, **kw):
        for p in (kw.pop("ps", None) or P_ARM):
            arms.append(dict(cell=f"{arm}_{level}_p{p}", arm=arm, level=str(level), N=N, K=10, p=p, **kw))
    for lz in (10.0, 20.0): add("A1", f"lamz{int(lz)}", lam_z=lz)
    add("A2", "maf0.01-0.05", maf=(0.01, 0.05))
    for g, nm in ((np.log(2), "gammaLog2"), (np.log(4), "gammaLog4")): add("A3", nm, gamma=g)
    add("A4", "prev0.01", design="casecontrol", prev=0.01, ps=[0.5, 0.2])
    for mlt in (0.5, 2.0, 5.0): add("A5", f"Nx{mlt}", n_mult=mlt)
    for fe in (0.2, 0.5, "flip"): add("A6", f"freqerr_{fe}", freq_err=fe)
    for h in ("matched", "severe"): add("A7", f"nhet_{h}", nhet=h)
    for inf in (0.8, 0.5): add("A8", f"info{inf}", info=inf)
    return core, arms


def seed_plan():
    core, arms = build_cells()
    root = np.random.SeedSequence(20260924)
    streams = root.spawn(1 + len(core) + len(arms))       # tuning stream spawned first
    tuning_streams = streams[0].spawn(len(core))
    plan = []
    for c, s in zip(core, tuning_streams):
        plan.append(("tuning", c, 200, s))
    for c, s in zip(core + arms, streams[1:]):
        plan.append(("test", c, 1000 if c["arm"] == "core" else 500, s))
    return plan


def one_replicate(setname, cell, rep, seq):
    rng = np.random.default_rng(seq)
    K = cell["K"]
    gcell = dict(cell); gcell["K"] = 3 * K                  # 3K candidates (Amendment 1)
    d = ss.generate(rng, gcell)
    out, fl = ss.analyse(d)
    lp = ss.pipeline_logp(out, fl["logp_spa"])
    sups = [("case_fraction", d["p_real"])]
    if cell.get("design") == "casecontrol":
        sups.append(("prevalence", cell["prev"]))
    rows = []
    for pl in ss.PIPELINES:
        idx, npass = ss.select(lp[pl], K)
        se = out[pl][1][idx]; f = d["f_rep"][idx]
        for supname, psup in sups:
            r = ss.classify(se, f, d["N_rep"], psup)
            r.update(set=setname, cell=cell["cell"], arm=cell["arm"], level=cell["level"],
                     p_target=cell["p"], K=K, rep=rep, pipeline=pl, true_scale=ss.TRUE_SCALE[pl],
                     p_supplied_as=supname, p_supplied=psup, p_real=d["p_real"], n_pass=npass,
                     conv_logit=bool(fl["P3"].all()), conv_firth=bool(fl["P4"].all()),
                     spa_ok_frac=float(fl["SPA"].mean()))
            rows.append(r)
    return rows


def run_chunk(task):
    setname, cell, start, stop, seqs = task
    fn = os.path.join(OUT, f"{setname}__{cell['cell']}__{start:04d}.parquet")
    if os.path.exists(fn):
        return fn, 0.0
    t = time.time()
    rows = []
    for i, s in zip(range(start, stop), seqs):
        rows += one_replicate(setname, cell, i, s)
    tmp = fn + ".tmp"
    pd.DataFrame(rows).to_parquet(tmp, index=False)
    os.replace(tmp, fn)
    return fn, time.time() - t


def tasks(only):
    T = []
    for setname, cell, reps, seq in seed_plan():
        if only != "all" and setname != only:
            continue
        kids = seq.spawn(reps)
        for a in range(0, reps, 50):
            T.append((setname, cell, a, min(a + 50, reps), kids[a:min(a + 50, reps)]))
    return T


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--only", default="all")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    T = tasks(a.only)
    T.sort(key=lambda t: -t[1]["K"] - (20 if t[1].get("design") == "casecontrol" else 0))  # longest first
    t0 = time.time(); done = 0
    with Pool(a.workers) as pool:
        for fn, sec in pool.imap_unordered(run_chunk, T):
            done += 1
            if done % 25 == 0 or done == len(T):
                print(f"{done}/{len(T)} chunks  elapsed {time.time()-t0:.0f}s", flush=True)
    print("complete")
