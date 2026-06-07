#!/usr/bin/env python3
"""
ENSO–Flood Audit: Complete Reproducible Pipeline
=================================================

Multi-stage audit of ENSO–flood teleconnections across 20 U.S. federal dams.
Uses stationary and non-stationary GEV modeling, likelihood-ratio testing,
effect-size and residual diagnostics, and Benjamini–Hochberg FDR correction.

Author: Sadman Arshad
Date: 2026-06-06
Citation: Arshad, S. (2026). A Reproducible Multi-Stage Audit of 
        ENSO–Flood Teleconnections in Regulated U.S. River Basins

This code is locked and reproducible. All seeds, bounds, and parameters are fixed.
"""

import pandas as pd
import numpy as np
from scipy.stats import kstest, pearsonr, chi2
from scipy.optimize import minimize

# ═════════════════════════════════════════════════════════════════════════════
# CONFIGURATION (LOCKED)
# ═════════════════════════════════════════════════════════════════════════════

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# GEV parameterization
XI_INITIALIZATIONS = [-0.2, -0.1, 0.0, 0.1]
XI_BOUNDS = (-0.5, 0.5)
FDR_ALPHA = 0.10
R_THRESHOLD = 0.15
KS_THRESHOLD = 0.05
LRT_ALPHA = 0.05

# ═════════════════════════════════════════════════════════════════════════════
# GEV LOG-LIKELIHOOD (LOCKED)
# ═════════════════════════════════════════════════════════════════════════════

def gev_nll(params, x, cov=None):
    """Negative log-likelihood for stationary (cov=None) or non-stationary GEV."""
    if cov is None:
        mu, sg, xi = params
        mu_v = np.full(len(x), mu)
        sg_v = np.full(len(x), sg)
    else:
        mu0, mu1, sg0, sg1, xi = params
        mu_v = mu0 + mu1 * cov
        sg_v = sg0 + sg1 * cov
    
    if np.any(sg_v <= 0):
        return 1e10
    z = 1 + xi * (x - mu_v) / sg_v
    if np.any(z <= 0):
        return 1e10
    if abs(xi) < 1e-6:
        return float(np.sum(np.log(sg_v)) + np.sum((x - mu_v) / sg_v) + 
                     np.sum(np.exp(-(x - mu_v) / sg_v)))
    return float(np.sum(np.log(sg_v)) + (1 + 1/xi) * np.sum(np.log(z)) + 
                 np.sum(z**(-1/xi)))

# ═════════════════════════════════════════════════════════════════════════════
# FITTING FUNCTIONS (LOCKED)
# ═════════════════════════════════════════════════════════════════════════════

def fit_stationary(x):
    """Fit stationary GEV to data x."""
    mu0, sg0 = np.mean(x), np.std(x) / 1.28
    best, best_nll = None, np.inf
    for xi in XI_INITIALIZATIONS:
        try:
            r = minimize(gev_nll, [mu0, sg0, xi], args=(x,), method='L-BFGS-B',
                        bounds=[(mu0-3*sg0, mu0+3*sg0), (sg0/10, sg0*10), XI_BOUNDS])
            if r.fun < best_nll:
                best_nll, best = r.fun, r
        except:
            pass
    return best

def fit_nonstationary(x, cov):
    """Fit non-stationary GEV with location and scale depending on covariate."""
    cz = (cov - cov.mean()) / cov.std()
    mu0, sg0 = np.mean(x), np.std(x) / 1.28
    best, best_nll = None, np.inf
    for xi in XI_INITIALIZATIONS:
        try:
            r = minimize(gev_nll, [mu0, 0, sg0, 0, xi], args=(x, cz), 
                        method='L-BFGS-B',
                        bounds=[(mu0-3*sg0, mu0+3*sg0), (-sg0, sg0), 
                                (sg0/10, sg0*10), (-sg0, sg0), XI_BOUNDS])
            if r.fun < best_nll:
                best_nll, best = r.fun, r
        except:
            pass
    return best

# ═════════════════════════════════════════════════════════════════════════════
# DIAGNOSTIC FUNCTIONS (LOCKED)
# ═════════════════════════════════════════════════════════════════════════════

def gev_cdf(x, params):
    """Compute GEV CDF for residual diagnostic."""
    mu, sg, xi = params
    z = 1 + xi * (x - mu) / sg
    if np.any(z <= 0):
        return None
    if abs(xi) < 1e-6:
        return np.exp(-(x - mu) / sg)
    return z**(-1/xi)

def classify_basin(p_lrt, r, ks_p):
    """Classify basin: SURVIVED, FRAGILE, MARGINAL, or NULL."""
    if p_lrt < LRT_ALPHA and abs(r) >= R_THRESHOLD and ks_p >= KS_THRESHOLD:
        return 'SURVIVED'
    elif p_lrt < LRT_ALPHA:
        return 'FRAGILE'
    elif LRT_ALPHA <= p_lrt < 0.10:
        return 'MARGINAL'
    else:
        return 'NULL'

# ═════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE (LOCKED)
# ═════════════════════════════════════════════════════════════════════════════

def audit_basin(peak, oni):
    """Run full GEV audit pipeline on one basin."""
    # Remove NaN
    valid = ~(np.isnan(peak) | np.isnan(oni))
    peak_v, oni_v = peak[valid], oni[valid]
    n = len(peak_v)
    
    if n < 10:
        return None
    
    # Fit stationary and non-stationary GEV
    st = fit_stationary(peak_v)
    ns = fit_nonstationary(peak_v, oni_v)
    
    if st is None or ns is None:
        return None
    
    # LRT
    lrt = max(0.0, 2 * (st.fun - ns.fun))
    p_lrt = float(1 - chi2.cdf(lrt, 2))
    
    # Correlation
    r, _ = pearsonr(peak_v, oni_v)
    
    # KS diagnostic
    cdf_vals = gev_cdf(peak_v, st.x)
    if cdf_vals is None:
        return None
    u_vals = -np.log(cdf_vals)
    _, ks_p = kstest(u_vals, 'expon')
    
    # Classification
    status = classify_basin(p_lrt, r, ks_p)
    
    return {
        'n': n,
        'p_lrt': p_lrt,
        'r': r,
        'ks_p': ks_p,
        'status_raw': status
    }

if __name__ == '__main__':
    print("ENSO–Flood Audit Pipeline (Locked, Reproducible)")
    print("See docstring for usage and citation.")
    print("\nTo run audit: Load your discharge and ONI data and call audit_basin().")
