#!/usr/bin/env python3
"""
Investigate the 10→5 projection (next in the 2n→n diagonal: 4→2, 6→3, 8→4, 10→5).
"""

import itertools
import math
import numpy as np
from scipy.linalg import det, orth, inv, sqrtm
from scipy.optimize import basinhopping
import sys

def normalize_tf(M, n, k):
    G = M @ M.T
    return np.real(inv(sqrtm(G))) @ M * math.sqrt(n / k)

def shadow_volume_orth(Q, n, k):
    return sum(abs(det(Q[list(idx)])) for idx in itertools.combinations(range(n), k))

np.random.seed(42)
n, k = 10, 5
print(f"Investigating {n}→{k} max-volume projection")
print(f"C({n},{k}) = {math.comb(n, k)} determinants per evaluation")

trials = 10
niter = 60
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
    sys.stdout.flush()

if best_P is None:
    print("Failed")
    sys.exit(1)

c = (best_P @ best_P.T)[0, 0]
gens = best_P.T / math.sqrt(c)
norms = np.linalg.norm(gens, axis=1)

print(f"\nVolume: {best_vol:.10f}")
print(f"Equilateral: {np.std(norms) < 1e-4} (std={np.std(norms):.6f})")
print(f"Tight frame: {np.allclose(best_P @ best_P.T, n/k * np.eye(k), atol=1e-4)}")

# Save
with open("results_10_5.txt", 'w') as f:
    f.write(f"10->5 max-volume projection results\n")
    f.write(f"Volume: {best_vol:.15f}\n")
    f.write(f"Equilateral: {np.std(norms) < 1e-4}\n")
    f.write(f"Norm std: {np.std(norms):.10f}\n\n")
    f.write(f"Generators ({n} vectors in R^{k}):\n")
    for i, (g, nm) in enumerate(zip(gens, norms)):
        coords = ", ".join(f"{g[j]:+.10f}" for j in range(k))
        f.write(f"  g{i:2d} = ({coords})  |g|={nm:.8f}\n")
print("Saved to results_10_5.txt")

# Symmetry (k=5, very slow for brute force — just check a subset)
print("\nComputing symmetry (limited search for k=5)...")
Ik = np.eye(k)
src_idx = None
for combo in itertools.combinations(range(n), k):
    if abs(det(gens[list(combo)])) > 1e-5:
        src_idx = list(combo)
        break

if src_idx:
    src = gens[src_idx]
    src_T_inv = np.linalg.inv(src.T)
    unique_R = []
    checked = 0
    # For k=5, P(10,5)*2^5 = 30240*32 ~ 1M — feasible but slow
    # Sample random permutations instead
    for _ in range(100000):
        perm = tuple(np.random.choice(n, k, replace=False))
        signs = tuple(np.random.choice([-1, 1], k))
        dst = np.array([signs[j] * gens[perm[j]] for j in range(k)])
        R = dst.T @ src_T_inv
        if not np.allclose(R @ R.T, Ik, atol=1e-3):
            continue
        if abs(abs(det(R)) - 1.0) > 1e-3:
            continue
        ok = all(any(np.allclose(R @ g, h, atol=1e-3) or np.allclose(R @ g, -h, atol=1e-3)
                     for h in gens) for g in gens)
        if ok and not any(np.allclose(R, S, atol=1e-3) for S in unique_R):
            unique_R.append(R)
            checked += 1

    print(f"Symmetry order (lower bound): {len(unique_R)}")
    orders = {}
    for R in unique_R:
        pw = Ik.copy()
        for o in range(1, 65):
            pw = pw @ R
            if np.allclose(pw, Ik, atol=1e-3):
                orders[o] = orders.get(o, 0) + 1
                break
    print(f"Element orders: {orders}")
    with open("results_10_5.txt", 'a') as f:
        f.write(f"\nSymmetry order (lower bound from random sampling): {len(unique_R)}\n")
        f.write(f"Element orders: {orders}\n")

# Inner products
G_inner = gens @ gens.T
ips = []
for i in range(n):
    for j in range(i+1, n):
        ips.append(round(G_inner[i, j], 5))
unique_ips = sorted(set(ips))
print(f"\nDistinct inner products: {len(unique_ips)}")
with open("results_10_5.txt", 'a') as f:
    f.write(f"\nDistinct inner products: {len(unique_ips)}\n")
    for ip in unique_ips:
        cnt = ips.count(ip)
        f.write(f"  {ip:+.6f}: {cnt} pairs\n")
        if cnt > 1:
            print(f"  {ip:+.6f}: {cnt} pairs  ***")

print("\nDone.")
