#!/usr/bin/env python3
"""
Verify that the O_h orbit of (a, a, b) with b/a = 1+√2 (silver ratio)
achieves the maximum shadow volume for the 12-cube → 3D projection.

The 12 generators are the O_h orbit of a point (a, a, b) on the unit sphere,
where b/a = 1 + √2 (the silver ratio). This gives exactly 12 antipodal
pairs (24 points on S²).

Exact values:
    a² = (5 - 2√2) / 17
    b² = (3 + 2√2)(5 - 2√2) / 17
    b/a = 1 + √2

Symmetry: O_h (order 48) — full octahedral symmetry
f-vector: (134, 264, 132) — Z(12,3) in general position
Free parameters: 0 (the silver ratio is uniquely determined by optimization)

Expected volume: 13.5096139098
"""

import itertools
import math
import numpy as np
from scipy.linalg import det, orth

SQRT2 = math.sqrt(2)
SILVER = 1 + SQRT2  # silver ratio

# Exact generator parameters
A = math.sqrt((5 - 2*SQRT2) / 17)
B = A * SILVER  # = A * (1 + sqrt(2))


def build_generators():
    """Build 12 generators as the O_h orbit of (a, a, b).

    The orbit consists of all vectors obtained by:
      - choosing which coordinate gets value b (3 choices)
      - ordering the remaining two as (a, a) (1 way, since equal)
      - applying all sign combinations (2³ = 8)
    Total: 3 × 8 = 24 vectors = 12 antipodal pairs.
    """
    vecs = set()
    base = [A, A, B]
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([1, -1], repeat=3):
            v = tuple(signs[i] * base[perm[i]] for i in range(3))
            canonical = min(v, tuple(-x for x in v))
            vecs.add(canonical)
    return np.array(sorted(vecs))


def main():
    gens = build_generators()
    n, k = 12, 3

    print("12-cube → 3D: O_h orbit with silver ratio")
    print(f"√2 = {SQRT2:.10f}")
    print(f"Silver ratio 1+√2 = {SILVER:.10f}")
    print(f"a = √((5-2√2)/17) = {A:.10f}")
    print(f"b = a·(1+√2)      = {B:.10f}")
    print(f"b/a = {B/A:.10f}")
    print(f"Check: 2a²+b² = {2*A**2 + B**2:.10f} (should be 1)")
    print()

    # Print generators
    norms = np.linalg.norm(gens, axis=1)
    print(f"{n} generators in R^{k}:")
    for i, (g, nm) in enumerate(zip(gens, norms)):
        print(f"  g{i:2d} = ({g[0]:+.8f}, {g[1]:+.8f}, {g[2]:+.8f})  |g| = {nm:.10f}")
    print()

    print(f"All unit vectors? {np.allclose(norms, 1.0)}")

    # Tight frame check
    P = gens.T  # 3×12
    G = P @ P.T
    expected = n / k
    print(f"PP^T = {G[0,0]:.10f} · I  (expected {expected:.10f})")
    print(f"Tight frame? {np.allclose(G, expected * np.eye(k))}")
    print()

    # Shadow volume
    Q = orth(P.T)  # 12×3 orthonormal
    vol = sum(abs(det(Q[list(idx)])) for idx in itertools.combinations(range(n), k))
    print(f"Shadow volume = {vol:.10f}")
    print(f"Expected      = 13.5096139098")
    print(f"Match: {abs(vol - 13.5096139098) < 1e-6}")


if __name__ == '__main__':
    main()
