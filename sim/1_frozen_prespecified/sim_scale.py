"""
Simulation of the three-signature effect-scale rule. Implements SIMULATION_PROTOCOL.md.

Every release is produced by fitting a regression to individual-level simulated data. No summary
statistic is generated from the signature formulas. All fits are vectorised across variants; a
mask matrix m (N x K, 0/1) restricts each variant's fit to its own analysed subset, which is how
variant-level sample-size heterogeneity (arm A7) is represented.
"""
import numpy as np
from scipy.special import expit, logit
from scipy.stats import norm

PIPELINES = ["P1", "P2", "P3", "P4", "P5", "P6", "P7"]
TRUE_SCALE = {"P1": "linear_probability", "P2": "standardised", "P3": "log_odds",
              "P4": "log_odds", "P5": "log_odds", "P6": "log_odds", "P7": "log_odds"}
SCALES = np.array(["linear_probability", "standardised", "log_odds"])


# ----------------------------------------------------------------------------- data generation
def _solve_alpha(eta, target, lo=-30.0, hi=30.0, it=80):
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if expit(mid + eta).mean() < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _studies_matched(rng, K):
    """Contributing-study counts matched to the CheckSumStats colorectal example
    (min 1, median 14, mean ~11.7, max 15)."""
    top = rng.choice([12, 13, 14, 15], size=K, p=[0.15, 0.15, 0.35, 0.35])
    low = rng.integers(1, 11, size=K)
    return np.where(rng.random(K) < 0.25, low, top)


def generate(rng, cell):
    """Returns a dict with genotype-for-analysis X_g (N,K), outcome y, covariate C, mask m,
    reported allele frequencies f_rep, reported N, realised case fraction, and truth fields."""
    N, K = cell["N"], cell["K"]
    lo, hi = cell.get("maf", (0.05, 0.50))
    f = rng.uniform(lo, hi, size=K)
    gamma = cell.get("gamma", 0.0)
    lam = cell.get("lam_z", 6.0)
    design = cell.get("design", "population")

    if design == "population":
        p = cell["p"]
        b = lam / np.sqrt(N * 2 * f * (1 - f) * p * (1 - p)) * rng.choice([-1, 1], size=K)
        G = rng.binomial(2, f, size=(N, K)).astype(float)
        C = rng.standard_normal(N)
        eta = G @ b + gamma * C
        a = _solve_alpha(eta, p)
        y = (rng.random(N) < expit(a + eta)).astype(float)
    else:  # case-control ascertainment from a population with prevalence cell["prev"]
        prev, cf = cell["prev"], cell["p"]
        n_case = int(round(cf * N)); n_ctrl = N - n_case
        b = lam / np.sqrt(N * 2 * f * (1 - f) * cf * (1 - cf)) * rng.choice([-1, 1], size=K)
        M = int(np.ceil(n_case / prev * 1.4))
        Gp = rng.binomial(2, f, size=(M, K)).astype(np.float32)
        Cp = rng.standard_normal(M).astype(np.float32)
        eta = Gp @ b.astype(np.float32) + np.float32(gamma) * Cp
        a = _solve_alpha(eta.astype(float), prev)
        yp = rng.random(M) < expit(a + eta)
        ic, io = np.flatnonzero(yp), np.flatnonzero(~yp)
        if len(ic) < n_case:
            raise RuntimeError("population too small for requested cases")
        idx = np.concatenate([rng.choice(ic, n_case, replace=False), rng.choice(io, n_ctrl, replace=False)])
        G = Gp[idx].astype(float); C = Cp[idx].astype(float)
        y = np.concatenate([np.ones(n_case), np.zeros(n_ctrl)])

    # imputation uncertainty: analysed dosage with squared correlation INFO to the true genotype
    info = cell.get("info", 1.0)
    if info < 1.0:
        e = rng.standard_normal((N, K)) * np.sqrt(2 * f * (1 - f))
        Xg = info * G + np.sqrt(info * (1 - info)) * e + (1 - info) * 2 * f
    else:
        Xg = G

    # variant-level sample-size heterogeneity
    m = np.ones((N, K))
    het = cell.get("nhet")
    if het is not None:
        if het == "matched":
            nj = np.round(N * _studies_matched(rng, K) / 15).astype(int)
        elif het == "severe":
            nj = np.where(np.arange(K) < K // 2, N // 4, N)
            rng.shuffle(nj)
        for j in range(K):
            keep = rng.permutation(N)[: nj[j]]
            m[:, j] = 0.0
            m[keep, j] = 1.0

    # reported allele frequency: from the analysed data, optionally perturbed (arm A6)
    f_rep = (m * Xg).sum(0) / (2 * m.sum(0))
    fe = cell.get("freq_err")
    if fe is not None:
        if fe == "flip":
            flip = rng.random(K) < 0.2
            f_rep = np.where(flip, 1 - f_rep, f_rep)
        else:
            f_rep = expit(logit(f_rep) + rng.normal(0, fe, size=K))

    N_rep = N * cell.get("n_mult", 1.0)
    p_real = y.mean()
    return dict(Xg=Xg, y=y, C=C, m=m, f_rep=f_rep, N_rep=N_rep, p_real=p_real, gamma=gamma)


# ----------------------------------------------------------------------------- estimators
def _design(Xg, C, use_cov):
    N, K = Xg.shape
    cols = [np.ones((N, K))]
    if use_cov:
        cols.append(np.repeat(C[:, None], K, 1))
    cols.append(Xg)
    return np.stack(cols, axis=2)  # (N, K, q), genotype is the last column


def fit_ols(y, X, m):
    q = X.shape[2]
    XtX = np.einsum("nk,nka,nkb->kab", m, X, X)
    Xty = np.einsum("nk,nka,n->ka", m, X, y)
    inv = np.linalg.inv(XtX)
    b = np.einsum("kab,kb->ka", inv, Xty)
    res = y[:, None] - np.einsum("nka,ka->nk", X, b)
    s2 = (m * res ** 2).sum(0) / (m.sum(0) - q)
    return b[:, -1], np.sqrt(s2 * inv[:, -1, -1])


def fit_logit(y, X, m, firth=False, maxit=60, tol=1e-9, maxstep=5.0):
    N, K, q = X.shape
    b = np.zeros((K, q))
    ybar = (m * y[:, None]).sum(0) / m.sum(0)
    b[:, 0] = logit(np.clip(ybar, 1e-6, 1 - 1e-6))
    conv = np.zeros(K, bool)
    for _ in range(maxit):
        mu = expit(np.einsum("nka,ka->nk", X, b))
        w = mu * (1 - mu)
        Iinv = np.linalg.inv(np.einsum("nk,nka,nkb->kab", m * w, X, X))
        r = y[:, None] - mu
        if firth:
            h = m * w * np.einsum("nka,kab,nkb->nk", X, Iinv, X)
            r = r + h * (0.5 - mu)
        U = np.einsum("nk,nka->ka", m * r, X)
        step = np.einsum("kab,kb->ka", Iinv, U)
        mx = np.abs(step).max(1)
        step = step / np.maximum(1.0, mx / maxstep)[:, None]
        b = b + step
        conv = mx < tol
        if conv.all():
            break
    mu = expit(np.einsum("nka,ka->nk", X, b))
    Iinv = np.linalg.inv(np.einsum("nk,nka,nkb->kab", m * mu * (1 - mu), X, X))
    return b[:, -1], np.sqrt(Iinv[:, -1, -1]), conv


def score_spa(y, Xg, C, m, use_cov, cutoff=2.0, maxit=60, tol=1e-10):
    """Covariate-adjusted score test with Lugannani-Rice saddlepoint P value for |z| > cutoff.
    Returns S, V, two-sided log P, and SPA convergence flags."""
    N, K = Xg.shape
    X0 = _design(np.zeros((N, K)), C, use_cov)[:, :, :-1]      # null design without genotype
    # null logistic fit (vectorised; differs by variant only through the mask)
    b0 = np.zeros((K, X0.shape[2]))
    ybar = (m * y[:, None]).sum(0) / m.sum(0)
    b0[:, 0] = logit(ybar)
    for _ in range(maxit):
        mu0 = expit(np.einsum("nka,ka->nk", X0, b0)); w0 = mu0 * (1 - mu0)
        I0 = np.einsum("nk,nka,nkb->kab", m * w0, X0, X0)
        st = np.linalg.solve(I0, np.einsum("nk,nka->ka", m * (y[:, None] - mu0), X0)[..., None])[..., 0]
        b0 = b0 + st
        if np.abs(st).max() < tol:
            break
    mu0 = expit(np.einsum("nka,ka->nk", X0, b0)); w0 = mu0 * (1 - mu0)
    A = np.einsum("nk,nka,nkb->kab", m * w0, X0, X0)
    B = np.einsum("nk,nka,nk->ka", m * w0, X0, Xg)
    gt = Xg - np.einsum("nka,ka->nk", X0, np.linalg.solve(A, B[..., None])[..., 0])
    gt = gt * m
    S = (gt * (y[:, None] - mu0)).sum(0)
    V = (gt ** 2 * w0).sum(0)
    z = S / np.sqrt(V)
    logp = np.log(2.0) + norm.logsf(np.abs(z))
    lmu = logit(mu0)
    lq0 = np.log1p(-mu0); lq1 = np.log(mu0)

    def cgf(t):
        x = t[None, :] * gt
        qq = expit(lmu + x)
        K0 = (np.logaddexp(lq0, lq1 + x) - x * mu0).sum(0)
        K1 = (gt * (qq - mu0)).sum(0)
        K2 = (gt ** 2 * qq * (1 - qq)).sum(0)
        return K0, K1, K2

    def tail(s):
        t = s / V
        ok = np.zeros(K, bool)
        for _ in range(maxit):
            _, K1, K2 = cgf(t)
            d = (K1 - s) / K2
            d = np.clip(d, -np.abs(t) - 1.0, np.abs(t) + 1.0)
            t = t - d
            ok = np.abs(d) < tol * (1 + np.abs(t))
            if ok.all():
                break
        K0, _, K2 = cgf(t)
        wv = np.sign(t) * np.sqrt(np.maximum(2 * (t * s - K0), 0.0))
        vv = t * np.sqrt(K2)
        with np.errstate(divide="ignore", invalid="ignore"):
            x = wv + np.log(vv / wv) / wv
        return x, ok

    use = np.abs(z) > cutoff
    spa_ok = np.ones(K, bool)
    if use.any():
        xu, oku = tail(np.abs(S))
        xl, okl = tail(-np.abs(S))
        lp = np.logaddexp(norm.logsf(xu), norm.logcdf(xl))
        good = use & oku & okl & np.isfinite(lp)
        logp = np.where(good, lp, logp)
        spa_ok = ~use | good
    return S, V, logp, spa_ok


def _z_from_logp(logp):
    return norm.isf(np.exp(logp - np.log(2.0)))


def analyse(d):
    """Apply P1-P7 to one simulated data set. Returns {pipeline: (beta, se)} and convergence info."""
    y, Xg, C, m = d["y"], d["Xg"], d["C"], d["m"]
    use_cov = d["gamma"] != 0.0
    X = _design(Xg, C, use_cov)
    out, flags = {}, {}
    out["P1"] = fit_ols(y, X, m)
    yz = (y - y.mean()) / y.std(ddof=1)
    out["P2"] = fit_ols(yz, X, m)
    b3, s3, c3 = fit_logit(y, X, m); out["P3"] = (b3, s3); flags["P3"] = c3
    b4, s4, c4 = fit_logit(y, X, m, firth=True); out["P4"] = (b4, s4); flags["P4"] = c4
    S, V, logp, ok = score_spa(y, Xg, C, m, use_cov); flags["SPA"] = ok
    zs = _z_from_logp(logp)
    out["P5"] = (b3, np.abs(b3) / zs)
    b6 = S / V
    out["P6"] = (b6, np.abs(b6) / zs)
    pr = d["p_real"]
    out["P7"] = (out["P1"][0] / (pr * (1 - pr)), out["P1"][1] / (pr * (1 - pr)))
    flags["logp_spa"] = logp
    return out, flags


def pipeline_logp(out, logp_spa):
    """Two-sided log P value of each pipeline's own test, used for instrument selection."""
    lp = {}
    for k in ["P1", "P2", "P3", "P4", "P7"]:
        b, s = out[k]
        lp[k] = np.log(2.0) + norm.logsf(np.abs(b / s))
    lp["P5"] = logp_spa
    lp["P6"] = logp_spa
    return lp


GWS_LOGP = np.log(5e-8)


def select(logp, K):
    """Indices of the first K candidates reaching genome-wide significance (candidates are in random order)."""
    idx = np.flatnonzero(np.isfinite(logp) & (logp < GWS_LOGP))
    return idx[:K], len(idx)


# ----------------------------------------------------------------------------- the rule under test
def classify(se, f_rep, N_rep, p_sup, tau=None, min_k=5):
    """Rules 1-6 of protocol section 3 for one release."""
    ok = np.isfinite(se) & (se > 0) & np.isfinite(f_rep) & (f_rep > 0) & (f_rep < 1)
    if ok.sum() < min_k:
        return dict(sigma_hat=np.nan, assigned="abstain_sufficiency", d_min=np.nan,
                    separation=np.nan, rel_sep=np.nan, d_lin=np.nan, d_std=np.nan, d_lo=np.nan)
    sh = np.median(se[ok] * np.sqrt(2 * N_rep * f_rep[ok] * (1 - f_rep[ok])))
    L = 0.5 * np.log(p_sup * (1 - p_sup))
    dist = np.abs(np.log(sh) - np.array([L, 0.0, -L]))
    o = np.argsort(dist)
    sep = dist[o[1]] - dist[o[0]]
    asg = SCALES[o[0]]
    if tau is not None and dist[o[0]] > tau:
        asg = "abstain_adequacy"
    return dict(sigma_hat=sh, assigned=asg, d_min=dist[o[0]], separation=sep, rel_sep=sep / abs(L),
                d_lin=dist[0], d_std=dist[1], d_lo=dist[2])
