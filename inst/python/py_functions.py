import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from numba import njit, prange

@njit(parallel=True)
def objective_numba(x, design, covar, group_start, group_end, case_offset):
    """Numba-compiled objective for conditional logistic WQS."""
    theta = x[0]
    n_exps = design.shape[1]
    total = 0.0
    # Parallel over matched sets
    for i in prange(group_start.shape[0]):
        s = group_start[i]
        e = group_end[i]
        co = case_offset[i]
        # linear predictor t_lin = design @ w + covar @ beta
        t_lin = np.empty(e - s, dtype=np.float64)
        for j in range(s, e):
            lin = 0.0
            # exposures
            for k in range(n_exps):
                lin += design[j, k] * x[1 + k]
            # covariates
            for k in range(covar.shape[1]):
                lin += covar[j, k] * x[1 + n_exps + k]
            t_lin[j - s] = lin
        # log-probabilities
        lps = theta * t_lin
        # stable log-sum-exp
        m = np.max(lps)
        lse = m + np.log(np.sum(np.exp(lps - m)))
        total += lps[co] - lse
    return -total

@njit
def gradient_numba(x, design, covar, group_start, group_end, case_offset):
    """Numba-compiled gradient of the objective (serial version to avoid reduction conflicts)."""
    theta = x[0]
    n_exps = design.shape[1]
    n_cov = covar.shape[1]
    g_theta = 0.0
    g_w = np.zeros(n_exps, dtype=np.float64)
    g_beta = np.zeros(n_cov, dtype=np.float64)
    # Serial over matched sets
    for i in range(group_start.shape[0]):
        s = group_start[i]
        e = group_end[i]
        co = case_offset[i]
        length = e - s
        # compute linear predictors and store
        t_lin = np.empty(length, dtype=np.float64)
        ds = np.empty((length, n_exps), dtype=np.float64)
        cs = np.empty((length, n_cov), dtype=np.float64)
        for j in range(s, e):
            lin = 0.0
            for k in range(n_exps):
                val = design[j, k]
                lin += val * x[1 + k]
                ds[j - s, k] = val
            for k in range(n_cov):
                val = covar[j, k]
                lin += val * x[1 + n_exps + k]
                cs[j - s, k] = val
            t_lin[j - s] = lin
        # softmax weights
        lps = theta * t_lin
        m = np.max(lps)
        exp_lps = np.exp(lps - m)
        sum_exp = np.sum(exp_lps)
        # accumulate gradient
        for j in range(length):
            p = exp_lps[j] / sum_exp
            # combine co and softmax in one reduction
            g_theta += (t_lin[j] if j == co else 0.0) - p * t_lin[j]
            for k in range(n_exps):
                g_w[k] += (ds[j, k] * theta if j == co else 0.0) - theta * p * ds[j, k]
            for k in range(n_cov):
                g_beta[k] += (cs[j, k] * theta if j == co else 0.0) - theta * p * cs[j, k]
    # assemble gradient vector
    grad = np.empty_like(x)
    grad[0] = -g_theta
    grad[1:1+n_exps] = -g_w
    grad[1+n_exps:] = -g_beta
    return grad


def pairwqs_noboot(wqsdata, col_vars, col_covars=None,
                          id_col="studyid", event_col="event", q=4):
    """
    Conditional logistic WQS (no bootstrap) using Numba for speed.
    Returns: dict with variable names, final weights, and wqs_beta.
    """
    import pandas as pd
    df = wqsdata.copy()
    # bin exposures into quantiles
    for col in col_vars:
        df[col] = pd.qcut(df[col], q=q, labels=False, duplicates='drop') + 1
    df['id'] = df[id_col]
    df['event'] = df[event_col].astype(int)

    # design & covariate matrices
    design = df[col_vars].to_numpy(dtype=np.float64).copy(order='C')
    if col_covars:
        covar = df[col_covars].to_numpy(dtype=np.float64).copy(order='C')
    else:
        covar = np.zeros((len(df), 0), dtype=np.float64)

    # sort by id
    ids = df['id'].values
    order = np.argsort(ids)
    design = design[order]
    covar = covar[order]
    event = df['event'].values[order]
    ids   = ids[order]

    # matched-set boundaries
    change = np.concatenate(([True], ids[1:] != ids[:-1]))
    group_start = np.nonzero(change)[0].astype(np.int64)
    group_end   = np.append(group_start[1:], len(ids)).astype(np.int64)
    # case offsets
    case_offset = np.empty(group_start.shape[0], dtype=np.int64)
    for i in range(group_start.shape[0]):
        s, e = group_start[i], group_end[i]
        grp = event[s:e]
        case_offset[i] = np.argmax(grp == 1)

    # initial parameters
    n_exps = design.shape[1]
    n_cov = covar.shape[1]
    x0 = np.concatenate(([0.2], np.full(n_exps, 1/n_exps), np.zeros(n_cov)))
    bounds = [(1e-6, None)] + [(1e-6, None)]*n_exps + [(None, None)]*n_cov
    cons = ({'type': 'eq', 'fun': lambda x: np.sum(x[1:1+n_exps]) - 1},)

    # optimize using JIT-compiled functions
    res = minimize(
        fun=lambda x: objective_numba(x, design, covar, group_start, group_end, case_offset),
        x0=x0,
        jac=lambda x: gradient_numba(x, design, covar, group_start, group_end, case_offset),
        method='SLSQP',
        bounds=bounds,
        constraints=cons,
        options={'maxiter': 500, 'ftol': 1e-6, 'disp': True}
    )

    return {"vars": col_vars,
            "final.weights": res.x[1:1+n_exps],
            "wqs_beta": res.x[0]}
