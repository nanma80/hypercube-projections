#!/usr/bin/env python3
"""
Compute the maximum-volume shadow of an n-cube projected to k dimensions.

Usage:
    python max_shadow_generic.py <n> <k> [--trials N] [--niter N]

Example:
    python max_shadow_generic.py 8 4
    python max_shadow_generic.py 10 3 --trials 30 --niter 200

This finds the k-dimensional subspace W ⊂ R^n that maximizes
Vol_k(proj_W([0,1]^n)), using basin-hopping global optimization
with tight-frame normalization.

Output:
    - Maximum shadow volume (normalized)
    - The k×n projection matrix P (tight frame, PP^T = (n/k)I)
    - The n generators in R^k (columns of P, normalized)
    - Generator norms and pairwise inner products
"""

import argparse
import itertools
import math
import sys
import numpy as np
from scipy.linalg import det, orth, inv, sqrtm
from scipy.optimize import minimize, basinhopping


def shadow_volume_orth(Q, n, k):
    """Shadow volume from n×k orthonormal matrix Q."""
    vol = 0.0
    for idx in itertools.combinations(range(n), k):
        vol += abs(det(Q[list(idx), :]))
    return vol


def shadow_volume(bases, n, k):
    """Shadow volume from a k×n projection matrix."""
    Q = orth(bases.T)  # n×k orthonormal
    return shadow_volume_orth(Q, n, k)


def normalize_to_tight_frame(M, n, k):
    """Normalize k×n matrix M to a tight frame: PP^T = (n/k)I_k."""
    G = M @ M.T  # k×k
    try:
        G_sqrt_inv = inv(sqrtm(G))
        P = np.real(G_sqrt_inv) @ M * math.sqrt(n / k)
    except Exception:
        P = M
    return P


def shadow_volume_tight(P, n, k):
    """Shadow volume from k×n tight frame matrix P."""
    c = (P @ P.T)[0, 0]
    gens = P.T / math.sqrt(c)  # n×k generators
    vol = 0.0
    for idx in itertools.combinations(range(n), k):
        vol += abs(det(gens[list(idx), :]))
    return vol


def maximize_shadow(n, k, n_trials=20, niter=150, verbose=True):
    """Find the maximum shadow volume for n-cube → k dimensions."""

    def neg_vol(params):
        M = params.reshape(k, n)
        P = normalize_to_tight_frame(M, n, k)
        return -shadow_volume_tight(P, n, k)

    best_vol = 0.0
    best_P = None

    for trial in range(n_trials):
        x0 = np.random.randn(k * n)
        try:
            res = basinhopping(neg_vol, x0, niter=niter,
                               minimizer_kwargs={'method': 'L-BFGS-B'},
                               disp=True, seed=trial)
            M = res.x.reshape(k, n)
            P = normalize_to_tight_frame(M, n, k)
            vol = shadow_volume_tight(P, n, k)
        except Exception:
            vol = 0.0
            P = None

        if vol > best_vol:
            best_vol = vol
            best_P = P
            if verbose:
                print(f"  Trial {trial+1}/{n_trials}: new best = {vol:.12f}")
        elif verbose and (trial + 1) % 5 == 0:
            print(f"  Trial {trial+1}/{n_trials}: vol = {vol:.12f}  (best = {best_vol:.12f})")

    return best_vol, best_P


def print_results(n, k, vol, P):
    """Print the optimization results."""
    c = (P @ P.T)[0, 0]
    gens = P.T / math.sqrt(c)  # n×k, rows are generators

    print()
    print("=" * 70)
    print(f"RESULT: {n}-cube → {k}D maximum shadow volume")
    print("=" * 70)
    print()

    # Volume
    print(f"Shadow volume: {vol:.15f}")
    print()

    # Tight frame check
    G = P @ P.T
    expected = n / k
    is_tight = np.allclose(G, expected * np.eye(k), atol=1e-6)
    print(f"Tight frame PP^T = {G[0,0]:.6f} · I_{k}  (expected {expected:.6f} · I_{k}, "
          f"{'✓' if is_tight else '✗'})")
    print()

    # Projection matrix
    print(f"Projection matrix P ({k}×{n}):")
    for i in range(k):
        row_str = "  [" + ", ".join(f"{P[i,j]:+.10f}" for j in range(n)) + "]"
        print(row_str)
    print()

    # Generators
    print(f"Generators (columns of P, normalized — {n} vectors in R^{k}):")
    norms = np.linalg.norm(gens, axis=1)
    for i in range(n):
        g = gens[i]
        coords = ", ".join(f"{g[j]:+.10f}" for j in range(k))
        print(f"  g{i:2d} = ({coords})  |g| = {norms[i]:.8f}")
    print()

    # Equilateral check
    equil = np.std(norms) < 1e-6
    print(f"Generator norms: min={norms.min():.8f}, max={norms.max():.8f}, "
          f"std={np.std(norms):.2e}")
    print(f"Equilateral (all norms equal): {'Yes' if equil else 'No'}")
    print()

    # Inner products
    print("Pairwise inner products (nonzero only):")
    G_inner = gens @ gens.T
    for i in range(n):
        for j in range(i+1, n):
            ip = G_inner[i, j]
            if abs(ip) > 1e-6:
                cos_angle = ip / (norms[i] * norms[j])
                print(f"  <g{i}, g{j}> = {ip:+.8f}  (cos = {cos_angle:+.8f})")
    print()

    # f-vector for the zonotope (general position formulas)
    print("Expected f-vector (general position zonotope Z(n,k)):")
    if k == 2:
        print(f"  Vertices: {2*n}, Edges: {2*n}, Faces: n/a (2D polygon)")
    elif k == 3:
        V = 2 * sum(math.comb(n-1, i) for i in range(k))
        E = n * 2 * sum(math.comb(n-2, i) for i in range(k-1))
        F = 2 * math.comb(n, 2)
        print(f"  Vertices: {V}, Edges: {E}, Faces: {F}")
        print(f"  Euler check: {V} - {E} + {F} = {V - E + F} (should be 2)")
    elif k == 4:
        V = 2 * sum(math.comb(n-1, i) for i in range(k))
        E = n * 2 * sum(math.comb(n-2, i) for i in range(k-1))
        F = math.comb(n, 2) * 2 * sum(math.comb(n-3, i) for i in range(k-2))
        C = 2 * math.comb(n, 3)
        print(f"  Vertices: {V}, Edges: {E}, 2-faces: {F}, 3-cells: {C}")
        print(f"  Euler check: {V} - {E} + {F} - {C} = {V - E + F - C} (should be 0)")


def main():
    parser = argparse.ArgumentParser(
        description="Compute maximum-volume shadow of n-cube projected to k dimensions.",
        epilog="Example: python max_shadow_generic.py 10 3 --trials 30")
    parser.add_argument("n", type=int, help="Dimension of the hypercube (high dimension)")
    parser.add_argument("k", type=int, help="Dimension of the projection (low dimension)")
    parser.add_argument("--trials", type=int, default=20,
                        help="Number of basin-hopping trials (default: 20)")
    parser.add_argument("--niter", type=int, default=150,
                        help="Basin-hopping iterations per trial (default: 150)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    args = parser.parse_args()

    n, k = args.n, args.k
    if k >= n:
        print(f"Error: k ({k}) must be less than n ({n})")
        sys.exit(1)
    if k < 1:
        print(f"Error: k must be at least 1")
        sys.exit(1)

    np.random.seed(args.seed)

    print(f"Computing max-volume shadow: {n}-cube → {k}D")
    print(f"  C({n},{k}) = {math.comb(n,k)} determinants per volume evaluation")
    print(f"  {args.trials} trials, {args.niter} iterations each")
    print()

    vol, P = maximize_shadow(n, k, n_trials=args.trials, niter=args.niter)

    if P is not None:
        print_results(n, k, vol, P)
    else:
        print("Optimization failed.")
        sys.exit(1)


if __name__ == '__main__':
    main()
