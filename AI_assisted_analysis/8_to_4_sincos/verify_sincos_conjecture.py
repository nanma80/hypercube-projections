#!/usr/bin/env python3
"""
Test the conjecture: the sin/cos parameterized matrix achieves the same
maximum shadow volume as a fully general 4×8 orthogonal projection.

The sin/cos family is a 2-parameter family within the 16-dimensional
Grassmannian Gr(4,8). We test whether the global optimum lies in this family.

Strategy:
  1. Run general (unconstrained) optimization many times with random starts.
  2. Run sin/cos optimization for comparison.
  3. Compare the best values found.
"""

import itertools
import math
import numpy as np
from scipy.linalg import orth, det
from scipy.optimize import minimize, basinhopping

# ── Shadow volume computation ────────────────────────────────────────────────

def shadow_volume_from_orth(Q):
    """Volume from an 8×4 orthonormal matrix Q (columns = basis of 4D subspace)."""
    vol = 0.0
    for idx in itertools.combinations(range(8), 4):
        vol += abs(det(Q[list(idx), :]))
    return vol

def shadow_volume_from_bases(bases):
    """Volume from a 4×8 matrix (rows span the subspace)."""
    Q = orth(bases.T)  # 8×4 orthonormal
    return shadow_volume_from_orth(Q)


# ── Method 1: General unconstrained optimization ────────────────────────────

def general_maximize(n_trials=30, niter_per_trial=200):
    """Optimize over ALL 4×8 matrices (full Grassmannian)."""
    
    def neg_vol_general(params):
        bases = params.reshape((4, 8))
        return -shadow_volume_from_bases(bases)
    
    best_vol = 0
    best_bases = None
    
    print(f"Running {n_trials} trials of general optimization (niter={niter_per_trial})...")
    for trial in range(n_trials):
        x0 = np.random.randn(32)  # 4×8 = 32 params
        try:
            res = basinhopping(neg_vol_general, x0, niter=niter_per_trial,
                              minimizer_kwargs={'method': 'L-BFGS-B'},
                              disp=False, seed=trial)
            vol = -res.fun
        except:
            vol = 0
        
        if vol > best_vol:
            best_vol = vol
            best_bases = res.x.reshape((4, 8))
            print(f"  Trial {trial+1}: NEW BEST = {vol:.12f}")
        elif trial % 10 == 0:
            print(f"  Trial {trial+1}: vol = {vol:.12f}, best = {best_vol:.12f}")
    
    return best_vol, best_bases


# ── Method 2: Sin/cos parameterized optimization ───────────────────────────

def get_bases_sincos(x, y):
    sx, cx = math.sin(x), math.cos(x)
    sy, cy = math.sin(y), math.cos(y)
    return np.array([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

def sincos_maximize(n_trials=50):
    """Optimize over the 2-parameter sin/cos family."""
    
    def neg_vol_sincos(params):
        return -shadow_volume_from_bases(get_bases_sincos(params[0], params[1]))
    
    best_vol = 0
    best_xy = None
    
    print(f"\nRunning {n_trials} trials of sin/cos optimization...")
    for trial in range(n_trials):
        x0 = np.random.rand(2) * 0.5
        res = basinhopping(neg_vol_sincos, x0, niter=100, disp=False, seed=trial)
        vol = -res.fun
        
        if vol > best_vol:
            best_vol = vol
            best_xy = res.x
    
    print(f"  Best sin/cos volume: {best_vol:.12f}")
    print(f"  Best angles: x={best_xy[0]:.15f}, y={best_xy[1]:.15f}")
    
    return best_vol, best_xy


# ── Method 3: Check if general optimum has the tight frame property ─────────

def check_tight_frame(bases):
    """Check if the projection has equal-length projections (tight frame)."""
    Q = orth(bases.T)  # 8×4
    diag = np.diag(Q @ Q.T)  # diagonal of QQ^T
    return diag, np.std(diag), np.allclose(diag, diag[0], atol=1e-6)


def check_sincos_structure(bases):
    """Check if a general optimal basis can be transformed to the sin/cos form."""
    Q = orth(bases.T)  # 8×4, orthonormal columns
    
    # The sin/cos structure implies specific relationships between columns.
    # Check: do the column norms of Q (as 4D vectors) match?
    # In the sin/cos form, all generators have norm 1/√2.
    norms = np.linalg.norm(Q, axis=1)
    
    # Check inner product structure
    G = Q @ Q.T  # 8×8 Gram-like matrix
    return G, norms


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    np.random.seed(42)
    
    print("=" * 70)
    print("CONJECTURE TEST: Does the sin/cos family achieve the global maximum?")
    print("=" * 70)
    print()
    print("Grassmannian Gr(4,8) has dimension 16.")
    print("Tight frame manifold has dimension ~9.")
    print("Sin/cos family has dimension 2.")
    print()
    
    # Method 2: Sin/cos (fast, do first)
    sincos_vol, sincos_xy = sincos_maximize(50)
    
    # Method 1: General optimization (slow, many trials)
    general_vol, general_bases = general_maximize(n_trials=30, niter_per_trial=200)
    
    # Compare
    print()
    print("=" * 70)
    print("COMPARISON")
    print("=" * 70)
    print(f"  Sin/cos max volume:  {sincos_vol:.15f}")
    print(f"  General max volume:  {general_vol:.15f}")
    print(f"  Difference:          {general_vol - sincos_vol:.2e}")
    print(f"  Relative difference: {(general_vol - sincos_vol)/sincos_vol:.2e}")
    
    if abs(general_vol - sincos_vol) < 1e-6:
        print(f"\n  ★ CONJECTURE SUPPORTED: Same maximum to within 1e-6.")
    elif general_vol > sincos_vol + 1e-6:
        print(f"\n  ✗ CONJECTURE REFUTED: General optimization found a higher volume!")
    else:
        print(f"\n  ? INCONCLUSIVE: Very close but not identical.")
    
    # Check tight frame property of general optimum
    if general_bases is not None:
        diag, std, is_tf = check_tight_frame(general_bases)
        print(f"\n  General optimum tight frame check:")
        print(f"    QQ^T diagonal: {np.round(diag, 6)}")
        print(f"    Std of diagonal: {std:.2e}")
        print(f"    Is tight frame: {is_tf}")
    
    # Also try: optimize on a larger scale with Nelder-Mead (more robust for noisy landscape)
    print("\n--- Additional verification with different optimizer ---")
    
    def neg_vol_general(params):
        bases = params.reshape((4, 8))
        return -shadow_volume_from_bases(bases)
    
    best_extra = 0
    for trial in range(20):
        x0 = np.random.randn(32) * 0.5
        res = minimize(neg_vol_general, x0, method='Nelder-Mead',
                      options={'maxiter': 50000, 'xatol': 1e-12, 'fatol': 1e-12})
        vol = -res.fun
        if vol > best_extra:
            best_extra = vol
    
    print(f"  Best Nelder-Mead general volume: {best_extra:.15f}")
    print(f"  vs sin/cos volume:               {sincos_vol:.15f}")
    print(f"  Difference: {best_extra - sincos_vol:.2e}")


if __name__ == '__main__':
    main()
