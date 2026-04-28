#!/usr/bin/env python3
"""
Compute the maximum-volume 2D shadow of the 4_21 polytope (E8 root system).

The 4_21 has 240 vertices: the E8 roots in R^8.
  - 112 integer roots: all permutations of (+/-1, +/-1, 0, 0, 0, 0, 0, 0)
  - 128 half-integer roots: (+/- 1/2)^8 with even number of minus signs

We search for a 2D subspace W of R^8 that maximizes Area(proj_W(4_21)).
"""
import argparse
import itertools
import numpy as np
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


def shadow_area(M):
    """M: 2x8 matrix. Project ROOTS via orthonormal basis, return hull area."""
    Q = orth(M.T)  # 8 x 2 orthonormal columns
    if Q.shape[1] < 2:
        return 0.0
    pts = ROOTS @ Q  # 240 x 2
    try:
        return ConvexHull(pts).volume  # in 2D, "volume" = area
    except Exception:
        return 0.0


def neg_area(params):
    return -shadow_area(params.reshape(2, 8))


def maximize_2d(n_trials=40, niter=100, seed=42, verbose=True):
    """Basin-hopping over 2x8 matrices to find max-area 2D projection."""
    rng = np.random.default_rng(seed)
    best_area = 0.0
    best_M = None
    for t in range(n_trials):
        x0 = rng.standard_normal(16)
        res = basinhopping(
            neg_area, x0, niter=niter,
            minimizer_kwargs={'method': 'L-BFGS-B',
                              'options': {'maxiter': 200, 'eps': 1e-4}},
            seed=int(rng.integers(1 << 30)),
            disp=False,
        )
        a = -res.fun
        if a > best_area:
            best_area = a
            best_M = res.x.reshape(2, 8)
            if verbose:
                print(f"  trial {t+1:3d}: new best = {a:.10f}")
        elif verbose and (t + 1) % 5 == 0:
            print(f"  trial {t+1:3d}: area = {a:.10f}  (best still {best_area:.10f})")
    return best_area, best_M


def get_coxeter_plane_basis():
    """The Coxeter plane projection (E8 Petrie polygon = 30-gon).

    Uses the two eigenvectors of the Coxeter element with eigenvalue
    exp(2*pi*i/30) — the standard 2D projection that shows E8's 30-fold
    symmetry.
    """
    b1 = []
    b2 = []
    for k in range(8):
        b1.append(np.cos(k * np.pi / 8 + np.pi / 16))
        b2.append(np.sin(k * np.pi / 8 + np.pi / 16))
    basis = np.array([b1, b2])
    return basis


def get_d8_petrie_basis():
    """D8 Petrie polygon basis (order-14 symmetry)."""
    b1 = []
    b2 = []
    for k in range(8):
        b1.append(np.cos(k * 2 * np.pi / 14))
        b2.append(np.sin(k * 2 * np.pi / 14))
    basis = np.array([b1, b2])
    return basis


def analyze_2d_projection(Q, label=""):
    """Analyze a 2D projection given orthonormal Q (8x2)."""
    pts = ROOTS @ Q  # 240 x 2
    hull = ConvexHull(pts)
    hv = pts[hull.vertices]
    n_hull = len(hull.vertices)

    print(f"\n{'=' * 60}")
    print(f"  {label}")
    print(f"{'=' * 60}")
    print(f"  Area: {hull.volume:.10f}")
    print(f"  Hull vertices: {n_hull}")
    print(f"  Hull edges: {len(hull.simplices)}")

    norms = np.linalg.norm(hv, axis=1)
    unique_norms = sorted(set(round(n, 4) for n in norms))
    print(f"  Hull vertex norms: {unique_norms}")
    equilateral = norms.std() < 1e-4
    print(f"  Equilateral: {'Yes' if equilateral else 'No'}")

    # Edge lengths
    edge_lengths = []
    for s in hull.simplices:
        d = np.linalg.norm(pts[s[0]] - pts[s[1]])
        edge_lengths.append(round(d, 4))
    unique_edges = sorted(set(edge_lengths))
    print(f"  Edge lengths: {unique_edges}")
    print(f"  Edge length counts: {[(l, edge_lengths.count(l)) for l in unique_edges]}")

    # Generator analysis
    gens = Q  # 8 x 2, row i = projection of e_i
    print(f"\n  8 generators (projections of standard basis vectors):")
    for i in range(8):
        g = gens[i]
        print(f"    g_{i} = ({g[0]:+.8f}, {g[1]:+.8f})  |g|={np.linalg.norm(g):.6f}")

    # Check coalescence
    print(f"\n  Generator coalescence check:")
    classes = []
    used = [False] * 8
    for i in range(8):
        if used[i]:
            continue
        cls = [i]
        used[i] = True
        for j in range(i + 1, 8):
            if not used[j]:
                if np.allclose(gens[i], gens[j], atol=1e-3):
                    cls.append(j)
                    used[j] = True
                elif np.allclose(gens[i], -gens[j], atol=1e-3):
                    cls.append(j)
                    used[j] = True
        classes.append(cls)
    for c in classes:
        g = gens[c[0]]
        print(f"    axes {c}: g = ({g[0]:+.6f}, {g[1]:+.6f})")
    if max(len(c) for c in classes) > 1:
        print(f"    => Coalescence detected!")
    else:
        print(f"    => No coalescence (all 8 generators distinct)")

    # Check for n-fold rotational symmetry
    print(f"\n  Symmetry check:")
    centroid = hv.mean(axis=0)
    angles = np.array([np.arctan2(v[1] - centroid[1], v[0] - centroid[0]) for v in hv])
    angles = np.sort(angles)
    diffs = np.diff(np.append(angles, angles[0] + 2 * np.pi))
    if np.std(diffs) < 1e-3:
        n_fold = n_hull
        print(f"    Regular {n_fold}-gon detected (equal angular spacing)")
    else:
        # Check for dihedral symmetry
        for n in range(2, 31):
            theta = 2 * np.pi / n
            R = np.array([[np.cos(theta), -np.sin(theta)],
                          [np.sin(theta), np.cos(theta)]])
            rotated = hv @ R.T
            # Check if rotated set matches original
            matched = True
            for p in rotated:
                if not any(np.allclose(p, q, atol=1e-3) for q in hv):
                    matched = False
                    break
            if matched:
                print(f"    C_{n} rotational symmetry detected")

    # Tight frame check
    QtQ = Q.T @ Q
    print(f"\n  Tight frame check: Q^T Q should be (8/2)*I_2 = 4*I_2")
    print(f"    Q^T Q = [[{QtQ[0,0]:.6f}, {QtQ[0,1]:.6f}], [{QtQ[1,0]:.6f}, {QtQ[1,1]:.6f}]]")

    return hull, hv, gens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=40)
    ap.add_argument('--niter', type=int, default=100)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    print(f"4_21 polytope: {len(ROOTS)} vertices in R^8\n")

    # ---- Baselines: known 2D projections ----
    print("BASELINES: Known 2D projections")

    # B8 Petrie
    b8_basis = np.array([
        [np.cos(k * np.pi / 8 + np.pi / 16) for k in range(8)],
        [np.sin(k * np.pi / 8 + np.pi / 16) for k in range(8)]
    ])
    Q_b8 = orth(b8_basis.T)
    analyze_2d_projection(Q_b8, label="B8 Petrie plane (16-fold)")

    # D8 Petrie
    d8_basis = get_d8_petrie_basis()
    Q_d8 = orth(d8_basis.T)
    analyze_2d_projection(Q_d8, label="D8 Petrie plane (14-fold)")

    # The phi-based 2D from e8.py (first two rows of get_bases)
    from math import sqrt
    phi = (1 + sqrt(5)) / 2
    e8_2d_basis = np.array([
        [1, phi, 0, -1, phi, 0, 0, 0],
        [phi, 0, 1, phi, 0, -1, 0, 0]
    ])
    Q_e8_2d = orth(e8_2d_basis.T)
    analyze_2d_projection(Q_e8_2d, label="E8.py phi-based 2D (first 2 rows of H4 basis)")

    # ---- Optimization ----
    print(f"\n{'=' * 60}")
    print(f"OPTIMIZING over 2D subspaces ({args.trials} trials, {args.niter} iters each)...")
    print(f"{'=' * 60}")
    best_area, best_M = maximize_2d(
        n_trials=args.trials, niter=args.niter, seed=args.seed)
    Q_opt = orth(best_M.T)  # 8 x 2

    hull_opt, hv_opt, gens_opt = analyze_2d_projection(
        Q_opt, label=f"MAX-AREA 2D PROJECTION (area={best_area:.10f})")

    # ---- Summary ----
    print(f"\n{'=' * 60}")
    print(f"  SUMMARY")
    print(f"{'=' * 60}")
    area_b8 = ConvexHull(ROOTS @ Q_b8).volume
    area_d8 = ConvexHull(ROOTS @ Q_d8).volume
    area_e8 = ConvexHull(ROOTS @ Q_e8_2d).volume
    print(f"  B8 Petrie area:     {area_b8:.10f}")
    print(f"  D8 Petrie area:     {area_d8:.10f}")
    print(f"  E8 phi-2D area:     {area_e8:.10f}")
    print(f"  Max-area:           {best_area:.10f}")
    print(f"  Ratio max/B8:       {best_area / area_b8:.6f}")
    print(f"  Ratio max/D8:       {best_area / area_d8:.6f}")


if __name__ == '__main__':
    main()
