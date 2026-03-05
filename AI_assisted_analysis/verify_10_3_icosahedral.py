#!/usr/bin/env python3
"""
Verify that the 10 three-fold axes of the icosahedron achieve the
maximum shadow volume for the 10-cube → 3D projection.

The 10 generators are the face-center directions of the icosahedron
(20 faces → 10 antipodal pairs → 10 directions). They form a tight
frame automatically: PP^T = (10/3)I_3.

Expected volume: 10.0840429737
"""

import itertools
import math
import numpy as np
from scipy.linalg import det

PHI = (1 + math.sqrt(5)) / 2  # golden ratio

# ── Exact generators ─────────────────────────────────────────────────────────
# N = sqrt((21 + 9*sqrt(5)) / 2)
#
# Type A (6 generators): cyclic permutations of (phi, 0, 2+sqrt(5)) / N
# Type B (4 generators): (±phi², ±phi², ±phi²) / N  (one per antipodal pair)

N = math.sqrt((21 + 9 * math.sqrt(5)) / 2)

generators = np.array([
    # Type A: cyclic permutations of (phi, 0, 2+sqrt(5)) / N
    [+PHI,      0,        2+math.sqrt(5)],   # g0
    [-PHI,      0,        2+math.sqrt(5)],   # g1
    [0,         2+math.sqrt(5),  +PHI],       # g4
    [0,         2+math.sqrt(5),  -PHI],       # g7
    [2+math.sqrt(5),  +PHI,      0],          # g8
    [2+math.sqrt(5),  -PHI,      0],          # g9
    # Type B: (±phi², ±phi², ±phi²) / N, one from each antipodal pair
    [+PHI**2,   +PHI**2,  +PHI**2],           # g2
    [-PHI**2,   +PHI**2,  +PHI**2],           # g3
    [+PHI**2,   +PHI**2,  -PHI**2],           # g5
    [-PHI**2,   +PHI**2,  -PHI**2],           # g6
]) / N


# ── Verify ───────────────────────────────────────────────────────────────────

def main():
    gens = generators
    n, k = 10, 3

    print("10-cube → 3D: Icosahedral three-fold axes")
    print(f"φ = (1+√5)/2 = {PHI:.10f}")
    print(f"N = √((21+9√5)/2) = {N:.10f}")
    print()

    # Print generators
    print(f"{n} generators in R^{k}:")
    norms = np.linalg.norm(gens, axis=1)
    for i, (g, nm) in enumerate(zip(gens, norms)):
        print(f"  g{i:2d} = ({g[0]:+.8f}, {g[1]:+.8f}, {g[2]:+.8f})  |g| = {nm:.10f}")
    print()

    # Unit vectors?
    print(f"All unit vectors? {np.allclose(norms, 1.0)}")

    # Tight frame check
    P = gens.T  # 3×10
    G = P @ P.T
    expected = n / k
    print(f"PP^T = {G[0,0]:.10f} · I  (expected {expected:.10f})")
    print(f"Tight frame? {np.allclose(G, expected * np.eye(k))}")
    print()

    # Shadow volume (using orthonormal basis, consistent with max_shadow_generic.py)
    from scipy.linalg import orth
    Q = orth(P.T)  # 10×3 orthonormal
    vol = sum(abs(det(Q[list(idx)])) for idx in itertools.combinations(range(n), k))
    print(f"Shadow volume = {vol:.10f}")
    print(f"Expected      = 10.0840429737")
    print(f"Match: {abs(vol - 10.0840429737) < 1e-6}")


if __name__ == '__main__':
    main()
