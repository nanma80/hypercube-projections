#!/usr/bin/env python3
"""
Investigate n→4 projections for arbitrary n.
Runs generic optimization, computes symmetry, saves results to file.

Usage: python investigate_n_4.py <n> [--trials N] [--niter N]
"""

import argparse
import itertools
import math
import numpy as np
from scipy.linalg import det, orth, inv, sqrtm
from scipy.optimize import basinhopping


def normalize_tf(M, n, k):
    G = M @ M.T
    return np.real(inv(sqrtm(G))) @ M * math.sqrt(n / k)


def shadow_volume_orth(Q, n, k):
    return sum(abs(det(Q[list(idx)])) for idx in itertools.combinations(range(n), k))


def optimize(n, k, trials, niter):
    best_vol = 0
    best_P = None

    for trial in range(trials):
        x0 = np.random.randn(k * n)
        def neg_vol(params):
            M = params.reshape(k, n)
            P = normalize_tf(M, n, k)
            Q = orth(P.T)
            return -shadow_volume_orth(Q, n, k)
        try:
            res = basinhopping(neg_vol, x0, niter=niter,
                               minimizer_kwargs={'method': 'L-BFGS-B'},
                               disp=False, seed=trial)
            M = res.x.reshape(k, n)
            P = normalize_tf(M, n, k)
            Q = orth(P.T)
            vol = shadow_volume_orth(Q, n, k)
        except Exception:
            vol = 0
            P = None

        if vol > best_vol:
            best_vol = vol
            best_P = P
            print(f"  Trial {trial+1}/{trials}: NEW BEST = {vol:.10f}")
        elif (trial + 1) % 5 == 0:
            print(f"  Trial {trial+1}/{trials}: vol = {vol:.10f}  (best = {best_vol:.10f})")

    return best_vol, best_P


def compute_symmetry(gens, n, k, tol=1e-3):
    Ik = np.eye(k)
    src_idx = None
    for combo in itertools.combinations(range(n), k):
        if abs(det(gens[list(combo)])) > 1e-5:
            src_idx = list(combo)
            break
    if src_idx is None:
        return -1, {}

    src = gens[src_idx]
    src_T_inv = np.linalg.inv(src.T)
    unique_R = []

    for perm in itertools.permutations(range(n), k):
        for signs in itertools.product([1, -1], repeat=k):
            dst = np.array([signs[j] * gens[perm[j]] for j in range(k)])
            R = dst.T @ src_T_inv
            if not np.allclose(R @ R.T, Ik, atol=tol):
                continue
            if abs(abs(det(R)) - 1.0) > tol:
                continue
            ok = all(any(np.allclose(R @ g, h, atol=tol) or np.allclose(R @ g, -h, atol=tol)
                         for h in gens) for g in gens)
            if ok and not any(np.allclose(R, S, atol=tol) for S in unique_R):
                unique_R.append(R)

    orders = {}
    for R in unique_R:
        pw = Ik.copy()
        for o in range(1, 65):
            pw = pw @ R
            if np.allclose(pw, Ik, atol=tol):
                orders[o] = orders.get(o, 0) + 1
                break

    return len(unique_R), orders


def get_generators(P, n, k):
    c = (P @ P.T)[0, 0]
    return P.T / math.sqrt(c)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int)
    parser.add_argument("--trials", type=int, default=12)
    parser.add_argument("--niter", type=int, default=80)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    n, k = args.n, 4
    np.random.seed(args.seed)

    print(f"Investigating {n}→{k} max-volume projection")
    print(f"C({n},{k}) = {math.comb(n, k)} determinants per evaluation")
    print(f"{args.trials} trials, {args.niter} iterations each")
    print()

    vol, P = optimize(n, k, args.trials, args.niter)

    if P is None:
        print("Optimization failed.")
        return

    gens = get_generators(P, n, k)
    norms = np.linalg.norm(gens, axis=1)
    equil = np.std(norms) < 1e-4

    print(f"\nVolume: {vol:.10f}")
    print(f"Equilateral: {equil} (std={np.std(norms):.6f})")
    print(f"Tight frame: {np.allclose(P @ P.T, n/k * np.eye(k), atol=1e-4)}")

    # Save generators
    outfile = f"results_{n}_4.txt"
    with open(outfile, 'w') as f:
        f.write(f"{n}->4 max-volume projection results\n")
        f.write(f"Volume: {vol:.15f}\n")
        f.write(f"Equilateral: {equil}\n")
        f.write(f"Norm std: {np.std(norms):.10f}\n\n")
        f.write(f"Generators ({n} vectors in R^4):\n")
        for i, (g, nm) in enumerate(zip(gens, norms)):
            f.write(f"  g{i:2d} = ({g[0]:+.10f}, {g[1]:+.10f}, {g[2]:+.10f}, {g[3]:+.10f})  |g|={nm:.8f}\n")
        f.write(f"\nProjection matrix P ({k}x{n}):\n")
        for i in range(k):
            f.write("  [" + ", ".join(f"{P[i,j]:+.10f}" for j in range(n)) + "]\n")
    print(f"\nGenerators saved to {outfile}")

    # Symmetry
    print(f"\nComputing symmetry...")
    sym_order, elem_orders = compute_symmetry(gens, n, k)
    print(f"Symmetry order: {sym_order}")
    print(f"Element orders: {elem_orders}")

    with open(outfile, 'a') as f:
        f.write(f"\nSymmetry order: {sym_order}\n")
        f.write(f"Element orders: {elem_orders}\n")

    # Inner products
    G_inner = gens @ gens.T
    ips = []
    for i in range(n):
        for j in range(i+1, n):
            ips.append(round(G_inner[i, j], 5))
    unique_ips = sorted(set(ips))
    print(f"\nDistinct inner products: {len(unique_ips)}")
    with open(outfile, 'a') as f:
        f.write(f"\nDistinct inner products: {len(unique_ips)}\n")
        for ip in unique_ips:
            cnt = ips.count(ip)
            f.write(f"  {ip:+.6f}: {cnt} pairs\n")
            if cnt > 1:
                print(f"  {ip:+.6f}: {cnt} pairs  ***")

    # Generator profiles
    profiles = {}
    for i, g in enumerate(gens):
        p = tuple(sorted(np.round(np.abs(g), 3)))
        if p not in profiles:
            profiles[p] = []
        profiles[p].append(i)
    print(f"\nDistinct |abs| profiles: {len(profiles)}")
    with open(outfile, 'a') as f:
        f.write(f"\nDistinct |abs| profiles: {len(profiles)}\n")
        for p in sorted(profiles.keys()):
            idxs = profiles[p]
            f.write(f"  {p}: {len(idxs)} generators\n")

    print(f"\nAll results saved to {outfile}")


if __name__ == '__main__':
    main()
