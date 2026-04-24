#!/usr/bin/env python3
"""
Deep investigation: WHY does the sin/cos family contain the global optimum,
and are there OTHER optima with potentially higher symmetry?

Part A: Structural analysis of the sin/cos family
Part B: Search for alternative optima via symmetry group actions
Part C: Orbit analysis under the full symmetry group of the problem
"""

import itertools
import math
import numpy as np
from scipy.linalg import det, orth, polar
from scipy.optimize import minimize, basinhopping
from collections import defaultdict

# ── Core helpers ─────────────────────────────────────────────────────────────

def shadow_volume_orth(Q):
    """Volume from 8×4 orthonormal Q."""
    vol = 0.0
    for idx in itertools.combinations(range(8), 4):
        vol += abs(det(Q[list(idx), :]))
    return vol

def shadow_volume_general(M):
    """Volume from a general 4×8 matrix."""
    Q = orth(M.T)
    return shadow_volume_orth(Q)

def get_sincos_basis(x, y):
    sx, cx = math.sin(x), math.cos(x)
    sy, cy = math.sin(y), math.cos(y)
    return np.array([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

def retract_to_stiefel(M):
    U, _, Vt = np.linalg.svd(M, full_matrices=False)
    return U @ Vt

# Known optimal angles
X_OPT = 0.073480735017031270149
Y_OPT = 0.203285503186765959386

# ══════════════════════════════════════════════════════════════════════════════
# PART A: Structural Analysis — WHY does the sin/cos family work?
# ══════════════════════════════════════════════════════════════════════════════

def part_a_structure():
    print("=" * 70)
    print("PART A: STRUCTURAL ANALYSIS OF THE SIN/COS FAMILY")
    print("=" * 70)

    P = get_sincos_basis(X_OPT, Y_OPT)

    # ── A1: Block structure ──
    print("\n--- A1: Block Structure ---")
    print("P has a 2×2 block structure in the row index:")
    print("  Rows 0-1 depend only on x (the angle in the (1,2) block)")
    print("  Rows 2-3 depend only on y (the angle in the (3,4) block)")
    print()

    # Verify: PP^T = 4I
    G = P @ P.T
    print(f"PP^T = \n{np.round(G, 10)}")
    print(f"PP^T = 4I? {np.allclose(G, 4*np.eye(4))}")

    # ── A2: Column structure — what pattern do the 8 columns have? ──
    print("\n--- A2: Column Analysis ---")
    cols = P.T  # 8×4
    print("Columns of P^T (= generators in R^4):")
    for i, c in enumerate(cols):
        print(f"  g{i}: [{c[0]:+.6f}, {c[1]:+.6f}, {c[2]:+.6f}, {c[3]:+.6f}]")

    # All generators have the same norm?
    norms = np.linalg.norm(cols, axis=1)
    print(f"\nGenerator norms: {np.round(norms, 10)}")
    print(f"All equal? {np.allclose(norms, norms[0])} (norm = {norms[0]:.10f})")
    print(f"Expected: sqrt(2) = {math.sqrt(2):.10f}")

    # ── A3: The sign pattern — relationship to Hadamard ──
    print("\n--- A3: Sign Pattern Analysis ---")
    # Extract the sign/magnitude structure
    # For x close to 0: sx ≈ x, cx ≈ 1, so the matrix is approximately:
    #   [x, x, 1, 1, -x, -x, 1, 1]
    #   [1, 1, x, x,  1,  1,-x,-x]
    #   [y,-y, Y,-Y, -Y,  Y, y,-y]
    #   [Y,-Y, y,-y,  y, -y,-Y, Y]
    # where Y = cos(y) ≈ 1.
    # The sign structure of P at x=y=π/4 would be:
    P_signs = get_sincos_basis(math.pi/4, math.pi/4)
    s = 1/math.sqrt(2)
    print(f"P at x=y=π/4 (normalized by 1/√2):")
    print(np.round(P_signs / s, 2))

    # Compare to Hadamard H8
    H2 = np.array([[1,1],[1,-1]])
    H4 = np.kron(H2, H2)
    H8 = np.kron(H4, H2)
    print(f"\nHadamard H8 (first 4 rows):")
    print(H8[:4])

    # ── A4: The key insight — P decomposes into two orthogonal 2D rotations ──
    print("\n--- A4: Decomposition into Independent Rotations ---")
    # Define rotation matrices
    Rx = np.array([[math.sin(X_OPT), math.cos(X_OPT)],
                   [math.cos(X_OPT), math.sin(X_OPT)]])  # not standard rotation
    print("Top-left 2×2 blocks of P:")
    print(f"  Cols(0,1): [{P[0,0]:.6f}, {P[0,1]:.6f}] / [{P[1,0]:.6f}, {P[1,1]:.6f}]")
    print(f"  Cols(2,3): [{P[0,2]:.6f}, {P[0,3]:.6f}] / [{P[1,2]:.6f}, {P[1,3]:.6f}]")
    print(f"  Cols(4,5): [{P[0,4]:.6f}, {P[0,5]:.6f}] / [{P[1,4]:.6f}, {P[1,5]:.6f}]")
    print(f"  Cols(6,7): [{P[0,6]:.6f}, {P[0,7]:.6f}] / [{P[1,6]:.6f}, {P[1,7]:.6f}]")

    # Rows 0-1 use the SAME 2×2 blocks with different sign patterns.
    # The 8 columns pair up as: (0,1), (2,3), (4,5), (6,7)
    # Each pair shares the same 2D projection in the (row0, row1) plane
    # and differs in the (row2, row3) plane.
    print("\nThe 8 generators split into 4 PAIRS by their (row0,row1) projection:")
    for pair_idx, (i, j) in enumerate([(0,1), (2,3), (4,5), (6,7)]):
        ci, cj = cols[i], cols[j]
        print(f"  Pair {pair_idx}: g{i} & g{j}")
        print(f"    (row0,row1): ({ci[0]:+.6f},{ci[1]:+.6f}) vs ({cj[0]:+.6f},{cj[1]:+.6f})")
        print(f"    (row2,row3): ({ci[2]:+.6f},{ci[3]:+.6f}) vs ({cj[2]:+.6f},{cj[3]:+.6f})")
        same_12 = np.allclose(ci[:2], cj[:2])
        diff_34 = np.allclose(ci[2:], -cj[2:]) or not np.allclose(ci[2:], cj[2:])
        print(f"    Same in (0,1)? {same_12}, Different in (2,3)? {diff_34}")

    # ── A5: Why tight frame + block structure constrains to 2 params ──
    print("\n--- A5: Why the Block Structure Forces 2 Parameters ---")
    print("""
The sin/cos family has three key structural properties:

1. BLOCK ORTHOGONALITY: P = [A; B] where A is 2×8 (depends on x) and
   B is 2×8 (depends on y), with A·B^T = 0 automatically.

2. PAIRED COLUMNS: The 8 columns pair as (0,1),(2,3),(4,5),(6,7).
   Within each pair, the top-2 components are identical and the
   bottom-2 components differ by a sign pattern.

3. TIGHT FRAME: PP^T = 4I is automatic for all x,y.

These three properties together constrain the family from 16 dimensions
(Gr(4,8)) down to just 2 dimensions (x,y). The question is whether
the volume-maximizing subspace NECESSARILY has these properties.
""")

    # ── A6: Check if the block structure is forced by optimality ──
    print("--- A6: Testing if General Optima Have Block Structure ---")
    # Run several general optimizations and check if results decompose into blocks
    np.random.seed(123)
    for trial in range(10):
        x0 = np.random.randn(32)
        res = basinhopping(lambda p: -shadow_volume_general(p.reshape(4,8)),
                           x0, niter=100,
                           minimizer_kwargs={'method': 'L-BFGS-B'},
                           disp=False, seed=trial)
        M = res.x.reshape(4, 8)
        Q = orth(M.T)  # 8×4 orthonormal

        # Check the Gram matrix structure
        G = Q @ Q.T  # 8×8
        # In the sin/cos family, G should have a specific block structure
        # related to the pairing of columns

        # Check if there's a 4×4 block decomposition
        # The subspace should decompose as V = V1 ⊕ V2 where V1, V2 are 2D
        # Check: does the 4×4 projection matrix Q Q^T have rank-2 blocks?

        # Actually, check if we can find a rotation R in O(4) such that
        # R·Q has the block form [A 0; 0 B] with A being 4×2 and B being 4×2
        # i.e., the subspace decomposes as a direct sum of two 2D subspaces.

        # Eigendecomposition of Q Q^T restricted to various 2D planes
        P_opt = Q.T  # 4×8
        # Check correlation between rows 0,1 and rows 2,3
        cross = P_opt[:2, :] @ P_opt[2:, :].T  # 2×2
        cross_norm = np.linalg.norm(cross)
        vol = -res.fun

        if trial == 0:
            print(f"  {'Trial':>5} {'Volume':>16} {'Cross-block norm':>18} {'Is block?':>10}")
        block = "YES" if cross_norm < 0.1 else "no"
        print(f"  {trial+1:5d} {vol:16.10f} {cross_norm:18.10f} {block:>10}")

    print("""
Note: The cross-block norm measures how much the optimal 4D subspace
decomposes into two orthogonal 2D subspaces. A value near 0 means
it naturally has the 2×2 block structure of the sin/cos family.

However, the general optimizer finds an ARBITRARY basis for the same
subspace, so we need to find the RIGHT rotation to reveal the blocks.
""")

    # Better approach: use SVD of the 8×4 projection matrix restricted to 
    # complementary sets of coordinates
    print("--- A6b: Proper Block Decomposition Check ---")
    print("Testing if optimal subspace decomposes as V = V1 ⊕ V2 (2D + 2D)")
    print()

    for trial in range(5):
        x0 = np.random.randn(32)
        res = basinhopping(lambda p: -shadow_volume_general(p.reshape(4,8)),
                           x0, niter=150,
                           minimizer_kwargs={'method': 'L-BFGS-B'},
                           disp=False, seed=trial+200)
        M = res.x.reshape(4, 8)
        Q = orth(M.T)  # 8×4 orthonormal
        vol = -res.fun

        if vol < 7.84:
            continue

        # The 4D subspace is spanned by Q's columns.
        # Check: can we find an orthogonal decomposition V = V1 ⊕ V2?
        # This is equivalent to: does the 4×4 Gram matrix of Q restricted to
        # various coordinate planes have a nice structure?

        # More precisely: find R ∈ O(4) such that (QR) has block structure
        # [A 0; 0 B] in some coordinate grouping.
        
        # The Gram matrix G = QQ^T is an 8×8 projection matrix.
        G = Q @ Q.T
        
        # Check all ways to partition {0,...,7} into two sets of 4
        # and see if the cross-block of G vanishes for any partition
        best_cross = float('inf')
        best_partition = None
        for S in itertools.combinations(range(8), 4):
            T = [j for j in range(8) if j not in S]
            cross_block = G[np.ix_(list(S), T)]
            cn = np.linalg.norm(cross_block, 'fro')
            if cn < best_cross:
                best_cross = cn
                best_partition = (S, tuple(T))

        print(f"  Trial {trial+1}: vol={vol:.10f}, best cross-block={best_cross:.6f}")
        print(f"    Partition: {best_partition[0]} | {best_partition[1]}")

    # ── A7: The deep reason — connection to representation theory ──
    print("\n--- A7: Representation-Theoretic Perspective ---")
    print("""
STRUCTURAL EXPLANATION:

The projection problem has a natural symmetry group: the hyperoctahedral
group B_8 = (Z/2)^8 ⋊ S_8 acts on R^8 by signed permutations of coordinates.
This preserves the hypercube [0,1]^8.

The sin/cos family P(x,y) has the following remarkable property:

  P(x,y) = R_x ⊗ S₁ + R_y ⊗ S₂

where R_x, R_y are 2×2 rotation-like matrices and S₁, S₂ are specific
sign-pattern matrices in R^{2×8}. This is a TENSOR DECOMPOSITION.

The fact that the optimal subspace decomposes as a tensor product of
two 2D structures is related to the fact that:

  (a) The tight frame condition forces equal-norm generators
  (b) The volume function is a sum of |det| terms, which rewards
      generators being "spread out" in 4D
  (c) Spreading generators optimally in 4D = R^2 × R^2 naturally
      leads to independent optimization in each R^2 factor

This reduces the 16-dimensional Grassmannian optimization to a
2-dimensional problem because the two R^2 factors decouple.
""")


# ══════════════════════════════════════════════════════════════════════════════
# PART B: Search for Alternative Optima
# ══════════════════════════════════════════════════════════════════════════════

def part_b_alternatives():
    print("\n" + "=" * 70)
    print("PART B: SEARCH FOR ALTERNATIVE OPTIMA")
    print("=" * 70)

    P_opt = get_sincos_basis(X_OPT, Y_OPT)
    Q_opt = orth(P_opt.T)  # 8×4
    vol_opt = shadow_volume_orth(Q_opt)
    print(f"\nReference volume: {vol_opt:.15f}")

    # ── B1: Action of signed permutations on the optimal subspace ──
    print("\n--- B1: Column Permutations (S_8 action) ---")
    print("Checking if permuting the 8 columns gives a DIFFERENT optimal subspace")
    print("with the same volume...")

    # For each permutation σ ∈ S_8, the permuted matrix P_σ has the same
    # column norms and inner product structure, but spans a DIFFERENT subspace
    # (in general). The volume may differ because |det| is not permutation-invariant.

    # We can't check all 8! = 40320 permutations efficiently,
    # but we can check the most natural ones.

    # The sin/cos structure pairs columns as (0,1),(2,3),(4,5),(6,7).
    # Permutations that preserve this pairing may give the same subspace.

    # Check: swapping pairs
    equiv_classes = defaultdict(list)

    # Check a representative set of permutations
    # First, the permutations that arise from the D_4×Z_2 symmetry group
    # (these should give the SAME subspace)
    # Then, others that might give DIFFERENT subspaces with the same volume.

    # Rather than enumerate S_8, let's check:
    # 1. Swap within pairs: (0↔1), (2↔3), (4↔5), (6↔7) — 2^4 = 16 sign changes
    # 2. Swap between pairs: e.g., (0,1)↔(2,3) — these exchange blocks
    # 3. Arbitrary permutations from a random sample

    unique_subspaces = []  # list of (Q, vol, description)

    def subspace_distance(Q1, Q2):
        """Principal angles between subspaces. 0 = same subspace."""
        M = Q1.T @ Q2  # 4×4
        sv = np.linalg.svd(M, compute_uv=False)
        sv = np.clip(sv, -1, 1)
        angles = np.arccos(sv)
        return np.max(angles)  # largest principal angle

    def is_new_subspace(Q_new, existing, tol=1e-6):
        for Q_ex, _, _ in existing:
            if subspace_distance(Q_new, Q_ex) < tol:
                return False
        return True

    # Start with the known optimum
    unique_subspaces.append((Q_opt, vol_opt, "original sin/cos"))

    # Check all 2^4 = 16 within-pair sign flips (negate some columns)
    print("\n  Checking within-pair column sign flips...")
    for signs in itertools.product([1, -1], repeat=4):
        P_new = P_opt.copy()
        for pair_idx, s in enumerate(signs):
            col = 2 * pair_idx
            P_new[:, col] *= s
            P_new[:, col+1] *= s
        Q_new = orth(P_new.T)
        vol_new = shadow_volume_orth(Q_new)
        if abs(vol_new - vol_opt) < 1e-6 and is_new_subspace(Q_new, unique_subspaces):
            unique_subspaces.append((Q_new, vol_new, f"pair signs {signs}"))

    print(f"  After pair sign flips: {len(unique_subspaces)} distinct subspaces")

    # Check pair permutations: all 4! = 24 ways to reorder the 4 pairs
    print("  Checking pair permutations...")
    pairs = [(0,1), (2,3), (4,5), (6,7)]
    for perm in itertools.permutations(range(4)):
        new_order = []
        for p in perm:
            new_order.extend(pairs[p])
        P_new = P_opt[:, new_order]
        Q_new = orth(P_new.T)
        vol_new = shadow_volume_orth(Q_new)
        if abs(vol_new - vol_opt) < 1e-6 and is_new_subspace(Q_new, unique_subspaces):
            unique_subspaces.append((Q_new, vol_new, f"pair perm {perm}"))

    print(f"  After pair permutations: {len(unique_subspaces)} distinct subspaces")

    # Check within-pair swaps: swap the two columns within each pair
    print("  Checking within-pair swaps...")
    for swaps in itertools.product([False, True], repeat=4):
        P_new = P_opt.copy()
        for pair_idx, do_swap in enumerate(swaps):
            if do_swap:
                c1, c2 = 2*pair_idx, 2*pair_idx+1
                P_new[:, [c1, c2]] = P_new[:, [c2, c1]]
        Q_new = orth(P_new.T)
        vol_new = shadow_volume_orth(Q_new)
        if abs(vol_new - vol_opt) < 1e-6 and is_new_subspace(Q_new, unique_subspaces):
            unique_subspaces.append((Q_new, vol_new, f"within-pair swaps {swaps}"))

    print(f"  After within-pair swaps: {len(unique_subspaces)} distinct subspaces")

    # Check combined: sign flips + pair perms + within-pair swaps
    # This is the full group preserving the pair structure
    print("  Checking combined pair symmetries...")
    count = 0
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product([1, -1], repeat=4):
            for swaps in itertools.product([False, True], repeat=4):
                new_order = []
                for p in perm:
                    c1, c2 = pairs[p]
                    if swaps[p]:
                        c1, c2 = c2, c1
                    new_order.extend([c1, c2])
                P_new = P_opt[:, new_order].copy()
                for pair_idx, s in enumerate(signs):
                    P_new[:, 2*pair_idx] *= s
                    P_new[:, 2*pair_idx+1] *= s
                Q_new = orth(P_new.T)
                vol_new = shadow_volume_orth(Q_new)
                if abs(vol_new - vol_opt) < 1e-6 and is_new_subspace(Q_new, unique_subspaces):
                    unique_subspaces.append((Q_new, vol_new, f"combined"))
                count += 1
    print(f"  Checked {count} transformations")
    print(f"  After combined: {len(unique_subspaces)} distinct subspaces")

    # ── B2: Random column permutations from full S_8 ──
    print("\n--- B2: Random S_8 Permutations + Sign Flips ---")
    np.random.seed(999)
    for trial in range(5000):
        perm = np.random.permutation(8)
        signs = np.random.choice([-1, 1], size=8)
        P_new = P_opt[:, perm] * signs[np.newaxis, :]
        Q_new = orth(P_new.T)
        vol_new = shadow_volume_orth(Q_new)
        if abs(vol_new - vol_opt) < 1e-6 and is_new_subspace(Q_new, unique_subspaces):
            unique_subspaces.append((Q_new, vol_new, f"random S8+signs trial {trial}"))

    print(f"  After 5000 random B_8 trials: {len(unique_subspaces)} distinct subspaces")

    # ── B3: O(4) left action — does rotating the 4D target give new optima? ──
    print("\n--- B3: Check if different sin/cos angles give same/different subspaces ---")
    # The critical point equations may have multiple solutions.
    # Check by solving the system from many starting points.
    def neg_vol_sincos(params):
        return -shadow_volume_general(get_sincos_basis(params[0], params[1]))

    sincos_optima = []
    for trial in range(200):
        x0 = np.random.rand(2) * math.pi - math.pi/2
        res = minimize(neg_vol_sincos, x0, method='Nelder-Mead',
                       options={'xatol': 1e-14, 'fatol': 1e-14, 'maxiter': 50000})
        vol = -res.fun
        if abs(vol - vol_opt) < 1e-6:
            xy = res.x
            P_new = get_sincos_basis(xy[0], xy[1])
            Q_new = orth(P_new.T)
            is_new = True
            for Q_ex, _, _ in sincos_optima:
                if subspace_distance(Q_new, Q_ex) < 1e-4:
                    is_new = False
                    break
            if is_new:
                sincos_optima.append((Q_new, vol, f"sin/cos x={xy[0]:.6f}, y={xy[1]:.6f}"))

    print(f"  Found {len(sincos_optima)} distinct sin/cos optima achieving max volume")
    for _, vol, desc in sincos_optima:
        print(f"    {desc}: vol={vol:.12f}")

    # Check if any of these are new subspaces not already found
    for Q_new, vol, desc in sincos_optima:
        if is_new_subspace(Q_new, unique_subspaces):
            unique_subspaces.append((Q_new, vol, desc))

    print(f"\n  Total distinct optimal subspaces found: {len(unique_subspaces)}")

    return unique_subspaces


# ══════════════════════════════════════════════════════════════════════════════
# PART C: Symmetry Analysis of All Optima
# ══════════════════════════════════════════════════════════════════════════════

def part_c_symmetry(unique_subspaces):
    print("\n" + "=" * 70)
    print("PART C: SYMMETRY ANALYSIS OF ALL OPTIMAL SUBSPACES")
    print("=" * 70)

    for idx, (Q, vol, desc) in enumerate(unique_subspaces[:10]):  # cap at 10
        print(f"\n--- Subspace {idx}: {desc} ---")
        print(f"  Volume: {vol:.12f}")

        # Compute generators in R^4
        P = Q.T  # 4×8
        gens = P.T  # 8×4

        # Find symmetries: orthogonal R ∈ O(4) such that R maps {±g_i} to {±g_j}
        symmetries = []
        I4 = np.eye(4)

        # Check all signed permutations of the 8 generators
        # A symmetry is R such that for every generator g_i,
        # R·g_i = ±g_{σ(i)} for some permutation σ and signs.

        # Fast approach: for each pair of 4-tuples of generators,
        # try to build R that maps one to the other.
        gen_indices = list(range(8))
        neg_gens = [-g for g in gens]

        # For a 4×4 orthogonal R, we need 4 independent constraints.
        # Pick 4 linearly independent generators and try all mappings.
        # Choose generators 0,1,2,3 (if they're independent)
        src = gens[:4]
        if abs(det(src)) < 1e-8:
            # Pick a different independent set
            for combo in itertools.combinations(range(8), 4):
                if abs(det(gens[list(combo)])) > 1e-8:
                    src = gens[list(combo)]
                    break

        src_inv = np.linalg.inv(src)

        # For each way to map {src generators} to {±generators}
        all_targets = list(range(8))
        target_options = []
        for i in range(8):
            target_options.append([(i, +1), (i, -1)])  # gen_i or -gen_i

        checked = 0
        for perm in itertools.permutations(range(8), 4):
            for signs in itertools.product([1, -1], repeat=4):
                dst = np.array([signs[k] * gens[perm[k]] for k in range(4)])
                R = dst @ src_inv
                # Check if R is orthogonal
                if not np.allclose(R @ R.T, I4, atol=1e-8):
                    continue
                # Check if R maps ALL 8 generators to ±generators
                is_sym = True
                for g in gens:
                    Rg = R @ g
                    found = False
                    for h in gens:
                        if np.allclose(Rg, h, atol=1e-8) or np.allclose(Rg, -h, atol=1e-8):
                            found = True
                            break
                    if not found:
                        is_sym = False
                        break
                if is_sym:
                    # Check if R is new
                    is_new_sym = True
                    for S in symmetries:
                        if np.allclose(R, S, atol=1e-8):
                            is_new_sym = False
                            break
                    if is_new_sym:
                        symmetries.append(R)
                checked += 1

        print(f"  Symmetry group order: {len(symmetries)}")

        # Classify element orders
        order_counts = {}
        for R in symmetries:
            power = I4.copy()
            for k in range(1, 33):
                power = power @ R
                if np.allclose(power, I4, atol=1e-8):
                    order_counts[k] = order_counts.get(k, 0) + 1
                    break

        print(f"  Element orders: {order_counts}")

        # Check if this is D_4 × Z_2 or something else
        if len(symmetries) == 16:
            print(f"  → D_4 × Z_2 (same as original)")
        elif len(symmetries) > 16:
            print(f"  → HIGHER SYMMETRY! ({len(symmetries)} > 16)")
        elif len(symmetries) < 16:
            print(f"  → Lower symmetry ({len(symmetries)} < 16)")

    # ── C2: Check pairwise relationships between optimal subspaces ──
    print("\n--- Pairwise Distances Between Optimal Subspaces ---")
    n = min(len(unique_subspaces), 10)
    for i in range(n):
        for j in range(i+1, n):
            d = subspace_distance_ext(unique_subspaces[i][0], unique_subspaces[j][0])
            print(f"  d(S{i}, S{j}) = {d:.6f} rad = {math.degrees(d):.2f}°")


def subspace_distance_ext(Q1, Q2):
    M = Q1.T @ Q2
    sv = np.linalg.svd(M, compute_uv=False)
    sv = np.clip(sv, -1, 1)
    return np.max(np.arccos(sv))


# ══════════════════════════════════════════════════════════════════════════════
# PART D: The Deeper Question — Can We PROVE the Block Decomposition?
# ══════════════════════════════════════════════════════════════════════════════

def part_d_proof_sketch():
    print("\n" + "=" * 70)
    print("PART D: EVIDENCE FOR BLOCK DECOMPOSITION THEOREM")
    print("=" * 70)

    print("""
CLAIM: The volume-maximizing 4D subspace for the 8-cube shadow
necessarily decomposes as V = V_1 ⊕ V_2 where V_1, V_2 ⊂ R^8 are
2-dimensional, and V_1, V_2 project onto complementary coordinate
pairs in a specific pattern.

EVIDENCE:

1. TIGHT FRAME NECESSITY (Ivanov 2018):
   Any local maximizer must be a tight frame: PP^T ∝ I_4.
   This means all 8 generators have equal norm in R^4.

2. VOLUME DECOMPOSITION:
   V = sum_{|S|=4} |det(P_S)| / (det PP^T)^{1/2}

   For a tight frame with PP^T = 4I, this simplifies to:
   V = (1/16) sum_{|S|=4} |det(P_S)|

   The 70 determinants can be grouped by which coordinates they use.

3. NUMERICAL EVIDENCE:
""")

    # Test if the optimal projection matrix can always be rotated into
    # the block form, regardless of how general optimization finds it.
    P_opt = get_sincos_basis(X_OPT, Y_OPT)
    Q_opt = orth(P_opt.T)  # 8×4

    # The Gram matrix G = Q Q^T is a rank-4 projection matrix in R^{8×8}
    G = Q_opt @ Q_opt.T
    print("   Gram matrix G = QQ^T (8×8 rank-4 projection):")
    print(f"   G =")
    for row in G:
        print(f"     [{', '.join(f'{v:+.6f}' for v in row)}]")

    # The block structure means G should have a 4+4 block decomposition
    # G = [G11  0 ; 0  G22] in some coordinate permutation
    # Check: is there a permutation of {0,...,7} such that G decomposes?

    print("\n   Checking for block-diagonal structure of G...")
    best_off_diag = float('inf')
    best_perm_partition = None

    # Check all (8 choose 4)/2 = 35 unordered partitions
    for S in itertools.combinations(range(8), 4):
        T = tuple(j for j in range(8) if j not in S)
        off_block = G[np.ix_(list(S), list(T))]
        off_norm = np.linalg.norm(off_block, 'fro')
        if off_norm < best_off_diag:
            best_off_diag = off_norm
            best_perm_partition = (S, T)

    print(f"   Best partition: {best_perm_partition[0]} | {best_perm_partition[1]}")
    print(f"   Off-diagonal block Frobenius norm: {best_off_diag:.10f}")
    if best_off_diag < 1e-6:
        print("   → G IS block-diagonal! The subspace decomposes as V_1 ⊕ V_2.")
    else:
        print("   → G is NOT block-diagonal in any coordinate partition.")
        print("   → The subspace does NOT decompose into coordinate-aligned 2D subspaces.")
        print("   → The block structure is in a ROTATED basis, not the standard basis.")

    # Check if the block structure appears after rotating within each pair
    print("\n   The sin/cos structure pairs coords as (0,1),(2,3),(4,5),(6,7).")
    print("   Checking 2×2 blocks of G for this pairing:")
    for i in range(4):
        for j in range(4):
            block = G[2*i:2*i+2, 2*j:2*j+2]
            print(f"   G[{2*i}:{2*i+2}, {2*j}:{2*j+2}] = {np.round(block, 6).tolist()}")

    # The actual block decomposition is in R^4 (target space), not R^8.
    # V = V_1 ⊕ V_2 means we can choose an ON basis for V such that
    # basis vectors 1,2 span V_1 and basis vectors 3,4 span V_2.
    print("\n   Checking R^4 block decomposition of general optimizers...")

    np.random.seed(42)
    for trial in range(5):
        x0 = np.random.randn(32)
        res = basinhopping(lambda p: -shadow_volume_general(p.reshape(4,8)),
                           x0, niter=200,
                           minimizer_kwargs={'method': 'L-BFGS-B'},
                           disp=False, seed=trial+500)
        vol = -res.fun
        if vol < 7.84:
            continue

        M = res.x.reshape(4, 8)
        Q = orth(M.T)  # 8×4

        # The 8 generators in R^4
        gens = Q  # 8×4

        # Check: can we find R ∈ O(4) that block-diagonalizes the generators?
        # i.e., find R such that R·gens has the pattern:
        # [a_i, b_i, 0, 0] for generators in group 1
        # [0, 0, c_i, d_i] for generators in group 2

        # Equivalent: find a 2D subspace W of R^4 such that 4 generators
        # lie in W and 4 lie in W^⊥.

        # Try all (8 choose 4) = 70 ways to split generators into two groups
        best_decomp = float('inf')
        best_split = None
        for S in itertools.combinations(range(8), 4):
            T = [j for j in range(8) if j not in S]
            # Generators in S should span a 2D subspace
            G_S = gens[list(S)]  # 4×4
            svs = np.linalg.svd(G_S, compute_uv=False)
            # If it's really 2D, last 2 singular values should be ~0
            residual = svs[2]**2 + svs[3]**2
            if residual < best_decomp:
                best_decomp = residual
                best_split = (S, tuple(T))

        print(f"   Trial {trial+1}: vol={vol:.10f}, best 2D split residual={best_decomp:.2e}, "
              f"partition={best_split[0]}|{best_split[1]}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    np.random.seed(42)
    part_a_structure()
    unique = part_b_alternatives()
    part_c_symmetry(unique)
    part_d_proof_sketch()

if __name__ == '__main__':
    main()
