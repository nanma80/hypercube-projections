#!/usr/bin/env python3
"""
Test the conjecture: the 15 two-fold axes of the icosahedron (edge-center
directions) achieve the maximum shadow volume for the 15-cube → 3D projection.

RESULT: CONJECTURE REFUTED.
  Icosahedral volume: 18.4970583145
  Numerical maximum:  19.0904463480
  Deficit: 3.11%

The pattern 6→3 (I_h, five-fold axes) and 10→3 (I_h, three-fold axes)
does NOT extend to 15→3 with two-fold axes. The icosahedral edge-center
directions are 3.1% below the true optimum.
"""

import itertools
import math
import numpy as np
from scipy.linalg import det, orth, inv, sqrtm
from scipy.optimize import basinhopping

PHI = (1 + math.sqrt(5)) / 2  # golden ratio


def build_icosahedral_15():
    """Build the 15 two-fold axis directions of the icosahedron."""
    # Icosahedron vertices: cyclic permutations of (0, ±1, ±φ)
    ico = []
    for cyc in [(0,1,2), (1,2,0), (2,0,1)]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = [0.0, 0.0, 0.0]
                v[cyc[1]] = s1
                v[cyc[2]] = s2 * PHI
                ico.append(v)
    ico = np.array(ico)

    # Find edges (vertex pairs at distance 2)
    edges = []
    for i in range(12):
        for j in range(i+1, 12):
            if abs(np.linalg.norm(ico[i] - ico[j]) - 2.0) < 0.01:
                edges.append((i, j))
    assert len(edges) == 30, f"Expected 30 edges, got {len(edges)}"

    # Edge midpoints as unit vectors
    midpoints = []
    for i, j in edges:
        m = (ico[i] + ico[j]) / 2
        midpoints.append(m / np.linalg.norm(m))
    midpoints = np.array(midpoints)

    # Pick one from each antipodal pair
    axes = []
    used = set()
    for i in range(len(midpoints)):
        if i in used:
            continue
        for j in range(i+1, len(midpoints)):
            if j not in used and np.allclose(midpoints[i], -midpoints[j], atol=1e-6):
                axes.append(midpoints[i])
                used.add(i)
                used.add(j)
                break
    assert len(axes) == 15, f"Expected 15 axes, got {len(axes)}"
    return np.array(axes), ico, edges


def shadow_volume(gens, n, k):
    """Shadow volume using orthonormal basis."""
    P = gens.T  # k×n
    Q = orth(P.T)  # n×k orthonormal
    return sum(abs(det(Q[list(idx)])) for idx in itertools.combinations(range(n), k))


def numerical_max(n=15, k=3, trials=15, niter=100):
    """Generic numerical optimization (no structural assumptions)."""
    def normalize_tf(M):
        G = M @ M.T
        return np.real(inv(sqrtm(G))) @ M * math.sqrt(n / k)

    best = 0
    for trial in range(trials):
        x0 = np.random.randn(k * n)
        def neg_vol(params):
            M = params.reshape(k, n)
            P = normalize_tf(M)
            Q = orth(P.T)
            return -sum(abs(det(Q[list(idx)])) for idx in itertools.combinations(range(n), k))
        try:
            res = basinhopping(neg_vol, x0, niter=niter,
                               minimizer_kwargs={'method': 'L-BFGS-B'},
                               disp=False, seed=trial)
            vol = -res.fun
            if vol > best:
                best = vol
                print(f"  Trial {trial+1}/{trials}: new best = {vol:.10f}")
            elif (trial+1) % 5 == 0:
                print(f"  Trial {trial+1}/{trials}: vol = {vol:.10f}  (best = {best:.10f})")
        except:
            pass
    return best


def main():
    gens, ico, edges = build_icosahedral_15()
    n, k = 15, 3

    print("15-cube → 3D: Icosahedral two-fold axes (edge centers)")
    print(f"φ = (1+√5)/2 = {PHI:.10f}")
    print()

    # Print generators
    norms = np.linalg.norm(gens, axis=1)
    print(f"{n} generators in R^{k}:")
    for i, (g, nm) in enumerate(zip(gens, norms)):
        print(f"  g{i:2d} = ({g[0]:+.8f}, {g[1]:+.8f}, {g[2]:+.8f})  |g| = {nm:.10f}")
    print()

    print(f"All unit vectors? {np.allclose(norms, 1.0)}")

    # Tight frame check
    P = gens.T  # 3×15
    G = P @ P.T
    expected = n / k
    print(f"PP^T = {np.round(G, 6)}")
    print(f"Expected: {expected:.6f} · I")
    print(f"Tight frame? {np.allclose(G, expected * np.eye(k))}")
    print()

    # Shadow volume from icosahedral axes
    vol_ico = shadow_volume(gens, n, k)
    print(f"Icosahedral shadow volume = {vol_ico:.10f}")
    print()

    # Numerical optimization for comparison
    print("Running numerical optimization (no structural assumptions)...")
    np.random.seed(42)
    vol_num = numerical_max(n, k, trials=15, niter=100)
    print()

    print("=" * 60)
    print(f"Icosahedral volume: {vol_ico:.10f}")
    print(f"Numerical maximum:  {vol_num:.10f}")
    diff = vol_num - vol_ico
    print(f"Difference:         {diff:+.2e}")
    if abs(diff) / max(vol_ico, vol_num) < 0.001:
        print("CONJECTURE SUPPORTED: icosahedral = numerical maximum")
    elif diff > 0:
        pct = diff / vol_num * 100
        print(f"CONJECTURE REFUTED: numerical is {pct:.2f}% higher")
    else:
        print("UNEXPECTED: icosahedral exceeds numerical (need more trials)")


if __name__ == '__main__':
    main()
