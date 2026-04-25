#!/usr/bin/env python3
"""
Compute the maximum-volume 3D shadow of the 4_21 polytope (E8 root system).

The 4_21 has 240 vertices: the E8 roots in R^8.
  - 112 integer roots: all permutations of (+/-1, +/-1, 0, 0, 0, 0, 0, 0)
  - 128 half-integer roots: (+/- 1/2)^8 with even number of minus signs

We search for a 3D subspace W of R^8 that maximizes Vol_3(proj_W(4_21)).

Compare with the previously-built 7-fold symmetric projection (D_7d), whose
volume we also compute as a baseline, and with the 4D max-volume result
(noted in e8.py: 4D max ~ 142.08, 4D 600-cell H4 ~ 129.44).
"""
import argparse
import itertools
import numpy as np
from scipy.linalg import orth
from scipy.optimize import basinhopping
from scipy.spatial import ConvexHull


def get_4_21_roots():
    roots = []
    # 112 integer roots: pick 2 of 8 positions, each +-1
    for i, j in itertools.combinations(range(8), 2):
        for si in (-1, 1):
            for sj in (-1, 1):
                v = [0.0] * 8
                v[i] = si
                v[j] = sj
                roots.append(v)
    # 128 half-integer roots: even number of minus signs
    for s in itertools.product((-1, 1), repeat=8):
        if int(np.prod(s)) == 1:
            roots.append([0.5 * x for x in s])
    return np.array(roots)  # 240 x 8


ROOTS = get_4_21_roots()


def shadow_volume(M):
    """M: 3x8 matrix. Project ROOTS via orthonormal basis from M, return hull volume."""
    Q = orth(M.T)  # 8 x 3 orthonormal columns
    if Q.shape[1] < 3:
        return 0.0
    pts = ROOTS @ Q  # 240 x 3
    try:
        return ConvexHull(pts).volume
    except Exception:
        return 0.0


def neg_vol(params):
    return -shadow_volume(params.reshape(3, 8))


def maximize(n_trials=20, niter=80, seed=42, verbose=True):
    rng = np.random.default_rng(seed)
    best_vol = 0.0
    best_M = None
    for t in range(n_trials):
        x0 = rng.standard_normal(24)
        res = basinhopping(
            neg_vol, x0, niter=niter,
            minimizer_kwargs={'method': 'L-BFGS-B',
                              'options': {'maxiter': 100, 'eps': 1e-4}},
            seed=int(rng.integers(1 << 30)),
            disp=False,
        )
        v = -res.fun
        if v > best_vol:
            best_vol = v
            best_M = res.x.reshape(3, 8)
            if verbose:
                print(f"  trial {t+1:3d}: new best = {v:.8f}")
        elif verbose:
            print(f"  trial {t+1:3d}: vol = {v:.8f}  (best = {best_vol:.8f})")
    return best_vol, best_M


# 7-fold projection from gen_421_to_3d.py: 8 generators in R^3 (physical coords)
def heptagonal_projection_basis():
    rho = 2 * np.cos(np.pi / 7)
    sigma = rho * rho - 1
    sigmaX2 = 2 * sigma
    skew = np.sin(3 * np.pi / 7)
    gens_vz = [
        [-4, -2, -2,  2,  0,  2, 1, 2, 1],
        [ 2, -2, -4,  0,  0, -2, 1, 2, 1],
        [-2, -2,  2, -4, -2, -2, 1, 2, 1],
        [ 2,  4,  0,  2, -2, -4, 1, 2, 1],
        [ 0,  2,  4, -2, -2,  2, 1, 2, 1],
        [ 2,  0,  2,  2,  4,  0, 1, 2, 1],
        [ 0,  0, -2,  0,  2,  4, 1, 2, 1],
        [ 0,  0,  0,  0,  0,  0, 5, 4, 1],
    ]
    G = []
    for v in gens_vz:
        x = v[0] + v[1] * rho + v[2] * sigma
        y = v[3] + v[4] * rho + v[5] * sigma
        z = v[6] + v[7] * rho + v[8] * sigma
        G.append([x + y / sigmaX2, y * skew, z])
    return np.array(G)  # 8 x 3 -- columns of P^T, rows are generators


def hept_volume():
    G = heptagonal_projection_basis()  # 8 x 3
    pts = ROOTS @ G  # 240 x 3
    return ConvexHull(pts).volume, G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=20)
    ap.add_argument('--niter', type=int, default=80)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    print(f"4_21: {len(ROOTS)} vertices in R^8")
    print()
    print("Baseline: 7-fold (D_7d) projection (heptagonal vZome model)")
    hv, G = hept_volume()
    print(f"  hull volume = {hv:.8f}")
    # The G is not orthonormal; we want to compare with the same input
    # normalization as the optimizer (orthonormal columns -> roots are
    # length sqrt(2) for ints / sqrt(2) for halves -> all length sqrt(2)).
    # But our 8-cube generators have ||g_k|| ~ 14.16; not unit-scale.
    # So we report two numbers:
    #   raw  : volume with original physical 8-cube generators
    #   norm : volume after projecting via the orthonormalized version
    Q = orth(G)  # 8 x 3 orthonormal cols spanning the same row-space as G^T... wait
    # G is 8x3 with rows = generators. The projection sends r in R^8 to
    # r @ G in R^3.  The image subspace is span(rows of G^T) = column space of G.
    # An orthonormal basis for that subspace is orth(G).
    Qn = orth(G)  # 8 x 3
    pts_norm = ROOTS @ Qn
    hv_norm = ConvexHull(pts_norm).volume
    print(f"  hull volume (orthonormalized basis) = {hv_norm:.8f}")
    print()

    print(f"Optimizing over 3D subspaces ({args.trials} trials, "
          f"{args.niter} iterations each)...")
    best_vol, best_M = maximize(n_trials=args.trials, niter=args.niter,
                                 seed=args.seed)
    print()
    print("=" * 70)
    print(f"MAX-VOLUME 4_21 -> 3D")
    print("=" * 70)
    print(f"Best volume: {best_vol:.10f}")
    print(f"Heptagonal baseline (orthonormalized): {hv_norm:.10f}")
    print(f"Ratio max / heptagonal = {best_vol / hv_norm:.6f}")
    print()
    Q = orth(best_M.T)  # 8 x 3
    print("Optimal 8x3 orthonormal projection matrix Q:")
    for i in range(8):
        print("  [" + ", ".join(f"{Q[i,j]:+.10f}" for j in range(3)) + "]")
    print()
    pts = ROOTS @ Q
    hull = ConvexHull(pts)
    print(f"Convex hull: {len(hull.vertices)} vertices, "
          f"{len(hull.simplices)} triangular facets")
    norms = np.linalg.norm(pts[hull.vertices], axis=1)
    print(f"Hull-vertex norms: min={norms.min():.6f}, max={norms.max():.6f}, "
          f"std={norms.std():.6f}")
    print(f"Equilateral hull vertices: "
          f"{'Yes' if norms.std() < 1e-4 else 'No'}")


if __name__ == '__main__':
    main()
