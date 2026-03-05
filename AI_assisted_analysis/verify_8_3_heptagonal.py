#!/usr/bin/env python3
"""
Verify the 8-cube → 3D max-shadow projection using the (7→2) + z-axis
decomposition.

The 8 generators decompose as:
  g_k = (r·cos(2πk/7), r·sin(2πk/7), h)   for k = 0,...,6
  g_7 = (0, 0, b)

where r = 4/√21 (fixed by tight frame), and h, b are determined by
the tight frame condition PP^T = (8/3)I and the volume optimization.

  r² = 16/21
  b² = 8/3 - 7h²
  h_opt ≈ 0.482556863960706

Symmetry: D_7 × Z_2 (order 28)
Expected volume: 7.0270375786
"""

import itertools
import math
import numpy as np
from scipy.linalg import det, orth
from scipy.optimize import minimize_scalar

# ── Parameterization ─────────────────────────────────────────────────────────
# r is fixed by tight frame. Only h is free.

R = 4 / math.sqrt(21)   # xy-radius, fixed
H_OPT = 0.482556863960706  # optimal z-height of 7 cone vectors


def build_generators(h=H_OPT):
    """Build the 8 generators for a given z-height h."""
    b_sq = 8/3 - 7 * h**2
    assert b_sq > 0, f"h={h} too large: b² = {b_sq} < 0"
    b = math.sqrt(b_sq)
    gens = []
    for k in range(7):
        phi = 2 * math.pi * k / 7
        gens.append([R * math.cos(phi), R * math.sin(phi), h])
    gens.append([0, 0, b])
    return np.array(gens), b


def shadow_volume(gens):
    """Shadow volume using orthonormal basis."""
    P = gens.T  # 3×8
    Q = orth(P.T)  # 8×3 orthonormal
    return sum(abs(det(Q[list(idx)])) for idx in itertools.combinations(range(8), 3))


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    # First verify by re-optimizing h
    def neg_vol(h):
        b_sq = 8/3 - 7 * h**2
        if b_sq <= 0:
            return 0
        gens, _ = build_generators(h)
        return -shadow_volume(gens)

    res = minimize_scalar(neg_vol, bounds=(0.01, math.sqrt(8/21) - 0.01),
                          method='bounded', options={'xatol': 1e-15})
    h_verified = res.x
    vol_verified = -res.fun

    # Build generators at the known optimum
    gens, b = build_generators(H_OPT)
    vol = shadow_volume(gens)
    a = math.sqrt(R**2 + H_OPT**2)
    theta = math.atan2(R, H_OPT)

    print("8-cube → 3D: (7→2) + z-axis decomposition")
    print(f"r = 4/√21 = {R:.10f} (fixed by tight frame)")
    print(f"h = {H_OPT:.15f} (optimal z-height)")
    print(f"b = √(8/3 - 7h²) = {b:.15f}")
    print(f"a = √(r²+h²) = {a:.15f} (cone generator norm)")
    print(f"θ = arctan(r/h) = {theta:.15f}")
    print()

    # Print generators
    norms = np.linalg.norm(gens, axis=1)
    print("8 generators in R^3:")
    for i, (g, nm) in enumerate(zip(gens, norms)):
        label = f"k={i}" if i < 7 else "axial"
        print(f"  g{i} = ({g[0]:+.10f}, {g[1]:+.10f}, {g[2]:+.10f})  |g| = {nm:.8f}  [{label}]")
    print()

    # Tight frame check
    P = gens.T
    G = P @ P.T
    print(f"PP^T = diag({G[0,0]:.10f}, {G[1,1]:.10f}, {G[2,2]:.10f})")
    print(f"Expected: (8/3)I = {8/3:.10f} · I")
    print(f"Tight frame? {np.allclose(G, 8/3 * np.eye(3))}")
    print()

    # Volume
    print(f"Shadow volume    = {vol:.10f}")
    print(f"Re-optimized h   = {h_verified:.15f}")
    print(f"Re-optimized vol = {vol_verified:.10f}")
    print(f"Expected         = 7.0270375786")
    print(f"Match: {abs(vol - 7.0270375786) < 1e-6}")


if __name__ == '__main__':
    main()
