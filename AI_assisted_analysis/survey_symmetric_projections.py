#!/usr/bin/env python3
"""
Systematic survey: For which n→k projections of the n-cube does the
volume-maximizing projection produce a HIGHLY SYMMETRIC zonotope?

Known results:
  n→2: ALWAYS the regular 2n-gon (D_n symmetry) — equally spaced vectors
  4→2: regular octagon, D_8 (order 16), silver ratio
  6→3: rhombic triacontahedron, I_h (order 120), golden ratio
  8→4: generic equilateral zonotope, D_4×Z_2 (order 16), degree-8 algebraic

We check:
  n→3 for n = 4, 5, 6, 7, 8, 9, 10, 12
  n→4 for n = 5, 6, 7, 8, 10, 12
"""

import itertools
import math
import numpy as np
from scipy.linalg import det, orth, sqrtm, inv
from scipy.optimize import minimize, basinhopping

PHI = (1 + math.sqrt(5)) / 2

# ── Core ─────────────────────────────────────────────────────────────────────

def shadow_volume_from_gens(gens, n, k):
    """Volume from n generators in R^k (n×k matrix, rows are generators)."""
    vol = 0.0
    for idx in itertools.combinations(range(n), k):
        vol += abs(det(gens[list(idx), :]))
    return vol

def normalize_to_tight_frame(M, n, k):
    """Given k×n matrix M, normalize to tight frame: P P^T = (n/k) I_k."""
    G = M @ M.T  # k×k
    try:
        G_sqrt_inv = inv(sqrtm(G))
        P = np.real(G_sqrt_inv) @ M * math.sqrt(n / k)
    except Exception:
        P = M
    return P

def generators_from_tight_frame(P):
    """Get generators from k×n tight frame P (PP^T = c I)."""
    k, n = P.shape
    c = (P @ P.T)[0, 0]
    return P.T / math.sqrt(c)  # n×k, each row is a generator

def shadow_volume_tight(P, n, k):
    """Volume from k×n tight frame matrix."""
    gens = generators_from_tight_frame(P)
    return shadow_volume_from_gens(gens, n, k)

def optimize_shadow(n, k, n_trials=30, niter=200):
    """Find the maximum shadow volume for projecting n-cube to k dimensions."""
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
                               disp=False, seed=trial)
            M = res.x.reshape(k, n)
            P = normalize_to_tight_frame(M, n, k)
            vol = shadow_volume_tight(P, n, k)
        except Exception:
            vol = 0.0
            P = None

        if vol > best_vol:
            best_vol = vol
            best_P = P

    return best_vol, best_P


def compute_symmetry_order(gens, n, k, tol=1e-4):
    """
    Compute symmetry group order of zonotope with given generators.
    A symmetry is R ∈ O(k) mapping the set {±g_i} to itself.
    
    gens is n×k (rows are generators, each a k-vector).
    R acts on column vectors: R @ g.T = result.T.
    """
    Ik = np.eye(k)

    # Find k independent generators
    src_idx = None
    for combo in itertools.combinations(range(n), k):
        sub = gens[list(combo)]
        if abs(det(sub)) > 1e-6:
            src_idx = list(combo)
            break
    if src_idx is None:
        return -1

    src = gens[src_idx]  # k×k, rows are source generators
    # We need R @ src[j,:].T = dst[j,:].T for each j
    # In matrix form: R @ src.T = dst.T, so R = dst.T @ inv(src.T)
    src_T_inv = np.linalg.inv(src.T)  # k×k

    unique_R = []
    for perm in itertools.permutations(range(n), k):
        for signs in itertools.product([1, -1], repeat=k):
            dst = np.array([signs[j] * gens[perm[j]] for j in range(k)])
            R = dst.T @ src_T_inv  # correct: R @ src.T = dst.T
            if not np.allclose(R @ R.T, Ik, atol=tol):
                continue
            if abs(abs(det(R)) - 1.0) > tol:
                continue
            ok = True
            for g in gens:
                Rg = R @ g  # R acts on row vector as R @ g (treating g as column)
                if not any(np.allclose(Rg, h, atol=tol) or np.allclose(Rg, -h, atol=tol)
                           for h in gens):
                    ok = False
                    break
            if ok and not any(np.allclose(R, S, atol=tol) for S in unique_R):
                unique_R.append(R)

    return len(unique_R)


def check_tight_frame(gens, n, k):
    """Check if generators form a tight frame (G G^T ∝ I)."""
    P = gens.T  # k×n
    G = P @ P.T  # k×k
    expected = G[0, 0]
    return np.allclose(G, expected * np.eye(k), atol=1e-4)


# ── Known solutions for validation ──────────────────────────────────────────

def icosahedral_6_3():
    """The 6→3 icosahedral solution (rhombic triacontahedron)."""
    # 6 face-diagonal directions of the icosahedron
    vecs = np.array([
        [0, 1, PHI],
        [0, 1, -PHI],
        [1, PHI, 0],
        [1, -PHI, 0],
        [PHI, 0, 1],
        [-PHI, 0, 1],
    ], dtype=float)
    # Normalize to tight frame: P P^T = 2 I_3
    P = vecs.T  # 3×6
    G = P @ P.T
    P = np.real(inv(sqrtm(G))) @ P * math.sqrt(2)
    return P

def sincos_8_4():
    """The 8→4 sin/cos solution."""
    x = 0.073480735017031270149
    y = 0.203285503186765959386
    sx, cx = math.sin(x), math.cos(x)
    sy, cy = math.sin(y), math.cos(y)
    return np.array([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])


def validate_known_cases():
    """Validate the symmetry computation on known solutions."""
    print("=" * 70)
    print("VALIDATION ON KNOWN SOLUTIONS")
    print("=" * 70)

    # 6→3 icosahedral
    P63 = icosahedral_6_3()
    gens63 = generators_from_tight_frame(P63)
    vol63 = shadow_volume_from_gens(gens63, 6, 3)
    tf63 = check_tight_frame(gens63, 6, 3)
    sym63 = compute_symmetry_order(gens63, 6, 3)
    print(f"\n6→3 icosahedral: vol={vol63:.10f}, tight={tf63}, sym={sym63}")
    print(f"  Expected: sym=120 (I_h)")

    # 8→4 sin/cos
    P84 = sincos_8_4()
    gens84 = generators_from_tight_frame(P84)
    vol84 = shadow_volume_from_gens(gens84, 8, 4)
    tf84 = check_tight_frame(gens84, 8, 4)
    sym84 = compute_symmetry_order(gens84, 8, 4)
    print(f"\n8→4 sin/cos: vol={vol84:.10f}, tight={tf84}, sym={sym84}")
    print(f"  Expected: sym=16 (D_4×Z_2)")


# ── Survey ───────────────────────────────────────────────────────────────────

def survey_projections(target_k, n_values):
    k = target_k
    max_sym = (2**k) * math.factorial(k)
    print(f"\n{'=' * 70}")
    print(f"SURVEY: n→{k} PROJECTIONS")
    print(f"Max possible symmetry: |B_{k}| = {max_sym}")
    print(f"{'=' * 70}")

    results = []

    for n in n_values:
        n_choose_k = math.comb(n, k)
        print(f"\n--- {n}→{k} ({n_choose_k} determinants) ---")

        # Scale effort to problem difficulty
        n_trials = max(15, 40 - 2*n)
        niter = max(80, 200 - 10*n)

        vol, P = optimize_shadow(n, k, n_trials=n_trials, niter=niter)

        if P is not None:
            gens = generators_from_tight_frame(P)
            tf = check_tight_frame(gens, n, k)
            norms = np.linalg.norm(gens, axis=1)
            equil = np.std(norms) < 1e-4

            print(f"  Volume: {vol:.10f}")
            print(f"  Tight frame: {tf}")
            print(f"  Equilateral: {equil} (norm={norms[0]:.6f})")

            if n <= 10:
                sym = compute_symmetry_order(gens, n, k)
                print(f"  Symmetry order: {sym}")
            else:
                sym = -1
                print(f"  Symmetry order: (skipped, n too large for brute force)")

            results.append((n, k, vol, sym, equil, tf))
        else:
            print(f"  FAILED")
            results.append((n, k, 0, -1, False, False))

    return results


def main():
    np.random.seed(42)

    print("SURVEY: Which n-cube → k-space max-volume projections")
    print("produce highly symmetric zonotopes?")
    print()

    # First validate
    validate_known_cases()

    # Then survey
    results_3d = survey_projections(3, [4, 5, 6, 7, 8, 9, 10])
    results_4d = survey_projections(4, [5, 6, 7, 8, 10, 12])

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY TABLE")
    print(f"{'=' * 70}")
    print(f"{'n→k':>6} {'Volume':>14} {'Sym':>6} {'Tight':>6} {'Equil':>6} {'Notes':>24}")
    print("-" * 70)

    # Known 2D results
    print(f"{'4→2':>6} {'4.828':>14} {'16':>6} {'Yes':>6} {'Yes':>6} {'Regular octagon (D_8)':>24}")
    print(f"{'6→2':>6} {'---':>14} {'24':>6} {'Yes':>6} {'Yes':>6} {'Regular 12-gon (D_12)':>24}")

    for n, k, vol, sym, equil, tf in results_3d + results_4d:
        equil_s = "Yes" if equil else "No"
        tf_s = "Yes" if tf else "No"
        sym_s = str(sym) if sym > 0 else "?"

        notes = ""
        if n == 6 and k == 3:
            notes = "Icosahedral (I_h)?"
        elif n == 8 and k == 4:
            notes = "D_4×Z_2"
        elif sym > 0 and sym >= 48 and k == 3:
            notes = "HIGH SYMMETRY!"
        elif sym > 0 and sym >= 48 and k == 4:
            notes = "HIGH SYMMETRY!"

        print(f"{f'{n}→{k}':>6} {vol:>14.6f} {sym_s:>6} {tf_s:>6} {equil_s:>6} {notes:>24}")


if __name__ == '__main__':
    main()
