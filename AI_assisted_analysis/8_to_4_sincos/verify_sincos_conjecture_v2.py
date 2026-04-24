#!/usr/bin/env python3
"""
Comprehensive verification of the sin/cos conjecture for 8D→4D max shadow.

Tests whether the 2-parameter sin/cos family achieves the global maximum
over the full 16-dimensional Grassmannian Gr(4,8).

Three independent verification methods:
  1. Basin-hopping with general 32-parameter optimization
  2. Nelder-Mead from many random starts
  3. Riemannian gradient descent on the Stiefel manifold
"""

import itertools
import math
import numpy as np
from scipy.linalg import det, orth, polar
from scipy.optimize import minimize, basinhopping

# ── Core computation ─────────────────────────────────────────────────────────

def shadow_volume_orth(Q):
    """Volume from 8×4 orthonormal Q. Sum of |det(4×4 submatrices)|."""
    vol = 0.0
    for idx in itertools.combinations(range(8), 4):
        vol += abs(det(Q[list(idx), :]))
    return vol

def shadow_volume_general(M):
    """Volume from a general 4×8 matrix (rows span subspace)."""
    Q = orth(M.T)  # 8×4 orthonormal
    return shadow_volume_orth(Q)

def get_sincos_basis(x, y):
    sx, cx = math.sin(x), math.cos(x)
    sy, cy = math.sin(y), math.cos(y)
    return np.array([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

def retract_to_stiefel(M):
    """Project a 4×8 matrix back to the Stiefel manifold (orthonormal rows)."""
    U, _, Vt = np.linalg.svd(M, full_matrices=False)
    return U @ Vt  # 4×8 with orthonormal rows

# ── Method 1: Basin-hopping (general 32-param) ──────────────────────────────

def method_basin_hopping(n_trials=40, niter=150):
    print("=" * 70)
    print("METHOD 1: Basin-hopping on general 4×8 matrices")
    print(f"  {n_trials} trials, {niter} iterations each")
    print("=" * 70)

    def neg_vol(params):
        M = params.reshape(4, 8)
        return -shadow_volume_general(M)

    best_vol = 0.0
    best_M = None

    for trial in range(n_trials):
        x0 = np.random.randn(32)
        try:
            res = basinhopping(neg_vol, x0, niter=niter,
                               minimizer_kwargs={'method': 'L-BFGS-B'},
                               disp=False, seed=trial)
            vol = -res.fun
        except Exception:
            vol = 0.0

        if vol > best_vol:
            best_vol = vol
            best_M = res.x.reshape(4, 8)
            print(f"  Trial {trial+1:3d}: NEW BEST = {vol:.12f}")
        elif (trial + 1) % 10 == 0:
            print(f"  Trial {trial+1:3d}: vol = {vol:.12f}  (best = {best_vol:.12f})")

    return best_vol, best_M

# ── Method 2: Nelder-Mead from random starts ────────────────────────────────

def method_nelder_mead(n_trials=30):
    print("\n" + "=" * 70)
    print("METHOD 2: Nelder-Mead from random starts")
    print(f"  {n_trials} trials")
    print("=" * 70)

    def neg_vol(params):
        M = params.reshape(4, 8)
        return -shadow_volume_general(M)

    best_vol = 0.0
    for trial in range(n_trials):
        x0 = np.random.randn(32) * 0.5
        res = minimize(neg_vol, x0, method='Nelder-Mead',
                       options={'maxiter': 50000, 'xatol': 1e-12, 'fatol': 1e-12})
        vol = -res.fun
        if vol > best_vol:
            best_vol = vol
            print(f"  Trial {trial+1:3d}: NEW BEST = {vol:.12f}")
        elif (trial + 1) % 10 == 0:
            print(f"  Trial {trial+1:3d}: vol = {vol:.12f}  (best = {best_vol:.12f})")

    return best_vol

# ── Method 3: Stiefel manifold optimization ──────────────────────────────────

def method_stiefel(n_trials=30, n_steps=500, lr=0.01):
    print("\n" + "=" * 70)
    print("METHOD 3: Projected gradient ascent on Stiefel manifold")
    print(f"  {n_trials} trials, {n_steps} steps each")
    print("=" * 70)

    best_vol = 0.0
    eps = 1e-6

    for trial in range(n_trials):
        # Random point on Stiefel manifold V_{4,8}
        M = np.random.randn(4, 8)
        M = retract_to_stiefel(M)

        for step in range(n_steps):
            vol0 = shadow_volume_general(M)
            grad = np.zeros((4, 8))
            for i in range(4):
                for j in range(8):
                    M[i, j] += eps
                    grad[i, j] = (shadow_volume_general(M) - vol0) / eps
                    M[i, j] -= eps
            M = M + lr * grad
            M = retract_to_stiefel(M)
            # Adaptive learning rate
            if step == n_steps // 2:
                lr *= 0.5

        vol = shadow_volume_general(M)
        if vol > best_vol:
            best_vol = vol
            print(f"  Trial {trial+1:3d}: NEW BEST = {vol:.12f}")
        elif (trial + 1) % 10 == 0:
            print(f"  Trial {trial+1:3d}: vol = {vol:.12f}  (best = {best_vol:.12f})")
        lr = 0.01  # reset

    return best_vol

# ── Sin/cos reference ────────────────────────────────────────────────────────

def sincos_optimum():
    print("\n" + "=" * 70)
    print("REFERENCE: Sin/cos family optimization")
    print("=" * 70)

    def neg_vol(params):
        M = get_sincos_basis(params[0], params[1])
        return -shadow_volume_general(M)

    best_vol = 0.0
    best_xy = None
    for trial in range(100):
        x0 = np.random.rand(2) * 0.5
        res = basinhopping(neg_vol, x0, niter=200, disp=False, seed=trial)
        vol = -res.fun
        if vol > best_vol:
            best_vol = vol
            best_xy = res.x

    print(f"  Best sin/cos volume:  {best_vol:.15f}")
    print(f"  Angles: x={best_xy[0]:.15f}, y={best_xy[1]:.15f}")
    return best_vol, best_xy

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    np.random.seed(42)

    print("CONJECTURE: The 2-parameter sin/cos family achieves the global")
    print("maximum shadow volume for 8D→4D projection of the hypercube.")
    print()

    # Reference value
    sc_vol, sc_xy = sincos_optimum()

    # Three independent methods
    bh_vol, bh_M = method_basin_hopping(n_trials=40, niter=150)
    nm_vol = method_nelder_mead(n_trials=30)
    st_vol = method_stiefel(n_trials=20, n_steps=300)

    # Summary
    print("\n" + "=" * 70)
    print("FINAL COMPARISON")
    print("=" * 70)
    print(f"  Sin/cos optimum:        {sc_vol:.15f}")
    print(f"  Basin-hopping best:     {bh_vol:.15f}  (diff: {bh_vol - sc_vol:+.2e})")
    print(f"  Nelder-Mead best:       {nm_vol:.15f}  (diff: {nm_vol - sc_vol:+.2e})")
    print(f"  Stiefel gradient best:  {st_vol:.15f}  (diff: {st_vol - sc_vol:+.2e})")

    general_best = max(bh_vol, nm_vol, st_vol)
    tol = 1e-6

    if abs(general_best - sc_vol) < tol:
        print(f"\n  ★ CONJECTURE SUPPORTED: All methods agree to within {tol:.0e}")
        print(f"    No general optimizer found a volume exceeding the sin/cos maximum.")
    elif general_best > sc_vol + tol:
        print(f"\n  ✗ CONJECTURE REFUTED: General optimization exceeded sin/cos by {general_best - sc_vol:.2e}")
    else:
        print(f"\n  ? INCONCLUSIVE: difference = {general_best - sc_vol:.2e}")

    # Check tight-frame property of best general solution
    if bh_M is not None:
        Q = orth(bh_M.T)
        diag = np.diag(Q @ Q.T)
        print(f"\n  Tight frame check (general best):")
        print(f"    QQ^T diagonal: {np.round(diag, 6)}")
        print(f"    Std: {np.std(diag):.2e} (0 = perfect tight frame)")

if __name__ == '__main__':
    main()
