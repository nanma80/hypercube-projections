#!/usr/bin/env python3
"""
Verify the maximum-volume 4D shadow of the 4_21 polytope (E8 root system).

Re-checks the results from the old e8.py (Python 2):
  - H4-symmetric (phi-based) projection -> volume ~129.44, 600 hull vertices
  - Numerical max -> volume ~142.08, 48 hull vertices, 144 edges

Also analyzes the optimal projection structure: generator coalescence,
symmetry group, etc.
"""
import argparse
import itertools
import numpy as np
from math import sqrt
from scipy.linalg import orth
from scipy.optimize import basinhopping
from scipy.spatial import ConvexHull


def get_4_21_roots():
    """Build the 240 E8 roots in R^8."""
    roots = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (-1, 1):
            for sj in (-1, 1):
                v = [0.0] * 8
                v[i] = si
                v[j] = sj
                roots.append(v)
    for s in itertools.product((-1, 1), repeat=8):
        if int(np.prod(s)) == 1:
            roots.append([0.5 * x for x in s])
    return np.array(roots)  # 240 x 8


ROOTS = get_4_21_roots()


def shadow_volume_4d(M):
    """M: 4x8 matrix. Project ROOTS via orthonormal basis, return hull volume."""
    Q = orth(M.T)  # 8 x 4 orthonormal columns
    if Q.shape[1] < 4:
        return 0.0
    pts = ROOTS @ Q  # 240 x 4
    try:
        return ConvexHull(pts).volume
    except Exception:
        return 0.0


def neg_vol_4d(params):
    return -shadow_volume_4d(params.reshape(4, 8))


def get_h4_basis():
    """The phi-based basis from e8.py that gives the 600-cell."""
    phi = (1 + sqrt(5)) / 2
    bases = [
        [1, phi, 0, -1, phi, 0, 0, 0],
        [phi, 0, 1, phi, 0, -1, 0, 0],
        [0, 1, phi, 0, -1, phi, 0, 0],
        [0, 0, 0, 0, 0, 0, phi + 1, phi - 1]
    ]
    return np.array(bases)


def get_edges_by_inner_product(vertices, tol=0.01):
    """Find edges: pairs of hull vertices with second-largest inner product."""
    n = len(vertices)
    if n < 2:
        return []
    inner_products = set()
    for i in range(min(n, 50)):
        for j in range(i + 1, min(n, 50)):
            ip = round(np.inner(vertices[i], vertices[j]), 4)
            inner_products.add(ip)
    inner_products = sorted(inner_products, reverse=True)
    if len(inner_products) < 2:
        return []
    target = inner_products[0]
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if abs(np.inner(vertices[i], vertices[j]) - target) < tol:
                edges.append((i, j))
    return edges


def maximize_4d(n_trials=30, niter=100, seed=42, verbose=True):
    """Basin-hopping over 4x8 matrices to find max-volume 4D projection."""
    rng = np.random.default_rng(seed)
    best_vol = 0.0
    best_M = None
    for t in range(n_trials):
        x0 = rng.standard_normal(32)
        res = basinhopping(
            neg_vol_4d, x0, niter=niter,
            minimizer_kwargs={'method': 'L-BFGS-B',
                              'options': {'maxiter': 200, 'eps': 1e-4}},
            seed=int(rng.integers(1 << 30)),
            disp=False,
        )
        v = -res.fun
        if v > best_vol:
            best_vol = v
            best_M = res.x.reshape(4, 8)
            if verbose:
                print(f"  trial {t+1:3d}: new best = {v:.8f}")
        elif verbose and (t + 1) % 5 == 0:
            print(f"  trial {t+1:3d}: vol = {v:.8f}  (best still {best_vol:.8f})")
    return best_vol, best_M


def analyze_projection(Q, label=""):
    """Analyze a 4D projection given orthonormal Q (8x4)."""
    pts = ROOTS @ Q  # 240 x 4
    hull = ConvexHull(pts)
    hv = pts[hull.vertices]

    print(f"\n{'=' * 60}")
    print(f"  {label}")
    print(f"{'=' * 60}")
    print(f"  Volume: {hull.volume:.10f}")
    print(f"  Hull vertices: {len(hull.vertices)}")
    print(f"  Hull simplices (3-faces): {len(hull.simplices)}")

    norms = np.linalg.norm(hv, axis=1)
    unique_norms = sorted(set(round(n, 4) for n in norms))
    print(f"  Hull vertex norms: {unique_norms}")
    equilateral = norms.std() < 1e-4
    print(f"  Equilateral: {'Yes' if equilateral else 'No'}")

    # Edge analysis via pairwise distances
    dists = []
    for i in range(len(hv)):
        for j in range(i + 1, len(hv)):
            dists.append(np.linalg.norm(hv[i] - hv[j]))
    dists = sorted(dists)
    min_dist = dists[0]

    # Edges = pairs at minimum distance
    edges = []
    for i in range(len(hv)):
        for j in range(i + 1, len(hv)):
            d = np.linalg.norm(hv[i] - hv[j])
            if d < min_dist * 1.05:
                edges.append((i, j))
    print(f"  Edges (nearest-neighbor): {len(edges)}")

    # Generator analysis (rows of Q = projections of standard basis vectors)
    print(f"\n  8 generators (rows of Q^T = columns of Q transposed):")
    gens = Q.T  # 8 x 4 (but Q is 8x4, so generators are *columns* of the 4x8 projection)
    # Actually Q is 8x4 orthonormal columns. The projection is ROOTS @ Q.
    # The i-th generator is Q[i, :] = the i-th row of Q.
    gens = Q  # 8 x 4, row i = projection of e_i
    for i in range(8):
        g = gens[i]
        print(f"    g_{i} = ({g[0]:+.8f}, {g[1]:+.8f}, {g[2]:+.8f}, {g[3]:+.8f})  |g|={np.linalg.norm(g):.6f}")

    # Check for generator coalescence
    print(f"\n  Generator coalescence check:")
    classes = []
    used = [False] * 8
    for i in range(8):
        if used[i]:
            continue
        cls = [i]
        used[i] = True
        for j in range(i + 1, 8):
            if not used[j] and np.allclose(gens[i], gens[j], atol=1e-3):
                cls.append(j)
                used[j] = True
        classes.append(cls)
    for c in classes:
        g = gens[c[0]]
        print(f"    axes {c}: g = ({g[0]:+.6f}, {g[1]:+.6f}, {g[2]:+.6f}, {g[3]:+.6f})")
    if max(len(c) for c in classes) > 1:
        print(f"    => Coalescence detected!")
    else:
        print(f"    => No coalescence (all 8 generators distinct)")

    # Check generator norms and Gram matrix
    gram = gens @ gens.T  # 8x8
    print(f"\n  Generator norms: {[round(np.linalg.norm(gens[i]), 6) for i in range(8)]}")
    print(f"  Tight frame check: Q^T Q should be (8/4)*I_4 = 2*I_4")
    QtQ = Q.T @ Q
    print(f"    Q^T Q diagonal: {[round(QtQ[i, i], 6) for i in range(4)]}")
    off_diag_max = max(abs(QtQ[i, j]) for i in range(4) for j in range(4) if i != j)
    print(f"    Q^T Q max off-diagonal: {off_diag_max:.2e}")

    return hull, hv, gens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=30)
    ap.add_argument('--niter', type=int, default=100)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    print(f"4_21 polytope: {len(ROOTS)} vertices in R^8\n")

    # ---- Baseline: H4 (phi-based) projection ----
    print("=" * 60)
    print("BASELINE: H4-symmetric (phi-based) projection from e8.py")
    print("=" * 60)
    H4_basis = get_h4_basis()  # 4 x 8
    Q_h4 = orth(H4_basis.T)  # 8 x 4 orthonormal
    pts_h4 = ROOTS @ Q_h4
    hull_h4 = ConvexHull(pts_h4)
    print(f"  Volume: {hull_h4.volume:.10f}")
    print(f"  Hull vertices: {len(hull_h4.vertices)}")
    hv_h4 = pts_h4[hull_h4.vertices]
    norms_h4 = np.linalg.norm(hv_h4, axis=1)
    print(f"  Hull vertex norms: min={norms_h4.min():.6f}, max={norms_h4.max():.6f}, std={norms_h4.std():.6f}")
    print(f"  Equilateral: {'Yes' if norms_h4.std() < 1e-4 else 'No'}")

    # Confirm: all 240 vertices on the hull? (for the 600-cell all should be)
    n_on_hull = len(hull_h4.vertices)
    print(f"  All 240 on hull: {'Yes' if n_on_hull == 240 else f'No ({n_on_hull}/240)'}")

    # Count edges via the standard inner product method from e8.py
    if n_on_hull <= 240:
        edges_h4 = get_edges_by_inner_product(hv_h4, tol=0.001)
        print(f"  Edges (inner product method): {len(edges_h4)}")

    # ---- Optimization: max-volume 4D projection ----
    print(f"\n{'=' * 60}")
    print(f"OPTIMIZING over 4D subspaces ({args.trials} trials, {args.niter} iters each)...")
    print(f"{'=' * 60}")
    best_vol, best_M = maximize_4d(
        n_trials=args.trials, niter=args.niter, seed=args.seed)
    Q_opt = orth(best_M.T)  # 8 x 4

    analyze_projection(Q_opt, label=f"MAX-VOLUME 4D PROJECTION (vol={best_vol:.8f})")

    # ---- Comparison ----
    print(f"\n{'=' * 60}")
    print(f"  COMPARISON")
    print(f"{'=' * 60}")
    print(f"  H4 (600-cell) volume:    {hull_h4.volume:.10f}")
    print(f"  Max-volume (optimized):  {best_vol:.10f}")
    print(f"  Ratio max/H4:            {best_vol / hull_h4.volume:.6f}")

    # ---- Also try the old e8.py optimal bases as a cross-check ----
    print(f"\n{'=' * 60}")
    print("CROSS-CHECK: e8.py's reported optimal bases")
    print(f"{'=' * 60}")
    old_optimal = np.array([
        [-0.20428321, -0.59858509,  0.06389811, -0.59506347,  0.24978869,
          0.05710992,  0.24540529, -0.34044246],
        [ 0.42419816, -0.46386136, -0.06555465, -0.18770753, -0.15416367,
         -0.46957864, -0.41506279,  0.38575818],
        [ 0.55479698, -0.11375784, -0.11333372, -0.08582672, -0.06081787,
          0.60916705,  0.40573804,  0.3458931 ],
        [ 0.27089508,  0.12824633,  0.76264977,  0.02096803,  0.55432917,
         -0.08949402,  0.00405637,  0.11308381]])
    Q_old = orth(old_optimal.T)  # 8 x 4
    pts_old = ROOTS @ Q_old
    hull_old = ConvexHull(pts_old)
    print(f"  Volume: {hull_old.volume:.10f}")
    print(f"  Hull vertices: {len(hull_old.vertices)}")
    analyze_projection(Q_old, label="OLD e8.py OPTIMAL BASES")


if __name__ == '__main__':
    main()
