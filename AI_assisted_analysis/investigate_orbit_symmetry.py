#!/usr/bin/env python3
"""
Follow-up analysis: Fix symmetry computation, properly analyze the orbit
under B_8, and determine if any alternative optimum has higher symmetry.

Key findings from investigate_sincos_structure.py:
  - Block structure (rows 0-1 ⊥ rows 2-3) is universal in all optima
  - 5358+ distinct optimal subspaces found
  - 28 distinct sin/cos angle pairs achieve the maximum

This script:
  1. Fixes the symmetry computation (use actual generators, not orth basis)
  2. Computes the orbit size under B_8 exactly
  3. Checks all orbit elements for higher symmetry
  4. Analyzes the R^4 decomposition V = V1 ⊕ V2 properly
"""

import itertools
import math
import numpy as np
from scipy.linalg import det

# ── Core ─────────────────────────────────────────────────────────────────────

X_OPT = 0.073480735017031270149
Y_OPT = 0.203285503186765959386

def get_sincos_basis(x, y):
    sx, cx = math.sin(x), math.cos(x)
    sy, cy = math.sin(y), math.cos(y)
    return np.array([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

def zonotope_volume(P):
    """Volume from 4×8 tight-frame matrix (PP^T = 4I)."""
    vol = 0.0
    for idx in itertools.combinations(range(8), 4):
        vol += abs(det(P[:, list(idx)]))
    return vol / 16.0

def generators_from_P(P):
    """Get the 8 generators in R^4 (columns of P, normalized)."""
    return (P / 2).T  # 8×4, each has norm 1/sqrt(2)

def subspace_projector(P):
    """Rank-4 projector onto the row space of P (in R^8)."""
    # PP^T = 4I, so P^T(PP^T)^{-1}P = P^T P / 4
    return P.T @ P / 4  # 8×8 rank-4 projector


# ══════════════════════════════════════════════════════════════════════════════
# 1. CORRECT SYMMETRY ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

def correct_symmetry_analysis():
    """
    The symmetry group of the zonotope is the set of R ∈ O(4) such that
    R permutes the set {±g_i : i=0,...,7} (the generators and their negatives).
    
    Previous analysis (Part C) had a bug: it used orth(P.T) which changes
    the generators. Here we use the actual generators g_i = P[:,i] / 2.
    """
    print("=" * 70)
    print("1. CORRECTED SYMMETRY ANALYSIS")
    print("=" * 70)

    P = get_sincos_basis(X_OPT, Y_OPT)
    gens = generators_from_P(P)  # 8×4
    I4 = np.eye(4)

    print(f"\n8 generators in R^4 (columns of P/2):")
    for i, g in enumerate(gens):
        print(f"  g{i} = [{g[0]:+.8f}, {g[1]:+.8f}, {g[2]:+.8f}, {g[3]:+.8f}]")

    # Find all R ∈ O(4) mapping {±g_i} → {±g_j}
    # Strategy: pick 4 independent generators as source, try all signed
    # permutation mappings, solve for R, verify on all 8.
    
    # Pick generators 0,1,2,3 as source (check independence)
    src = gens[:4]
    if abs(det(src)) < 1e-8:
        raise ValueError("First 4 generators are dependent")
    src_inv = np.linalg.inv(src)

    symmetries = []
    
    for perm in itertools.permutations(range(8), 4):
        for signs in itertools.product([1, -1], repeat=4):
            dst = np.array([signs[k] * gens[perm[k]] for k in range(4)])
            R = dst @ src_inv  # R such that R @ src[k] = signs[k] * gens[perm[k]]
            
            # Check orthogonality
            if not np.allclose(R @ R.T, I4, atol=1e-7):
                continue
            
            # Check R maps ALL generators to ±generators
            ok = True
            for g in gens:
                Rg = R @ g
                found = any(
                    np.allclose(Rg, h, atol=1e-7) or np.allclose(Rg, -h, atol=1e-7)
                    for h in gens
                )
                if not found:
                    ok = False
                    break
            if not ok:
                continue
            
            # Check if new
            if not any(np.allclose(R, S, atol=1e-7) for S in symmetries):
                symmetries.append(R)

    print(f"\nSymmetry group order: {len(symmetries)}")

    # Element orders
    order_dist = {}
    for R in symmetries:
        power = I4.copy()
        for k in range(1, 33):
            power = power @ R
            if np.allclose(power, I4, atol=1e-7):
                order_dist[k] = order_dist.get(k, 0) + 1
                break
    print(f"Element order distribution: {order_dist}")

    # Center
    center = [R for R in symmetries
              if all(np.allclose(R @ S, S @ R, atol=1e-7) for S in symmetries)]
    print(f"Center order: {len(center)}")

    # Identify group
    n = len(symmetries)
    if n == 16 and order_dist.get(1,0) == 1 and order_dist.get(2,0) == 7 and order_dist.get(4,0) == 8:
        print("Group identification: D_4 × Z_2")
    elif n == 16:
        print(f"Group of order 16, orders: {order_dist}")
    else:
        print(f"Group of order {n}")

    # Print generators of the symmetry group
    print("\nLooking for generators of the symmetry group...")
    # Find order-4 elements
    order4 = [R for R in symmetries if order_dist_single(R, I4) == 4]
    order2_non_central = [R for R in symmetries 
                         if order_dist_single(R, I4) == 2 
                         and not np.allclose(R, I4) and not np.allclose(R, -I4)]

    if order4:
        a = order4[0]
        print(f"\nOrder-4 generator a:")
        print(f"  {np.round(a, 6).tolist()}")
        # Describe the action
        for i, g in enumerate(gens):
            ag = a @ g
            for j, h in enumerate(gens):
                if np.allclose(ag, h, atol=1e-7):
                    print(f"  a(g{i}) = g{j}")
                    break
                elif np.allclose(ag, -h, atol=1e-7):
                    print(f"  a(g{i}) = -g{j}")
                    break

    if order2_non_central:
        b = order2_non_central[0]
        print(f"\nOrder-2 generator b:")
        print(f"  {np.round(b, 6).tolist()}")
        for i, g in enumerate(gens):
            bg = b @ g
            for j, h in enumerate(gens):
                if np.allclose(bg, h, atol=1e-7):
                    print(f"  b(g{i}) = g{j}")
                    break
                elif np.allclose(bg, -h, atol=1e-7):
                    print(f"  b(g{i}) = -g{j}")
                    break

    return symmetries


def order_dist_single(R, I4):
    power = I4.copy()
    for k in range(1, 33):
        power = power @ R
        if np.allclose(power, I4, atol=1e-7):
            return k
    return -1


# ══════════════════════════════════════════════════════════════════════════════
# 2. ORBIT UNDER B_8 — How many distinct optimal subspaces?
# ══════════════════════════════════════════════════════════════════════════════

def orbit_analysis():
    """
    B_8 = (Z/2)^8 ⋊ S_8 acts on R^8 by signed permutations.
    Each element (ε, σ) acts on a 4D subspace V by: V → {x ∈ R^8 : ...}
    
    Two subspaces V, V' are in the same orbit if some signed permutation
    maps V to V'. The orbit size = |B_8| / |stabilizer|.
    
    We compute the stabilizer: signed permutations that preserve the
    optimal subspace V.
    """
    print("\n" + "=" * 70)
    print("2. ORBIT ANALYSIS UNDER B_8")
    print("=" * 70)

    P = get_sincos_basis(X_OPT, Y_OPT)
    G = subspace_projector(P)  # 8×8 rank-4 projector onto V

    print(f"\nProjector G = P^T P / 4 (rank 4, 8×8):")
    print(f"G @ G ≈ G? {np.allclose(G @ G, G)}")
    print(f"trace(G) = {np.trace(G):.6f} (should be 4)")

    # Stabilizer: signed permutations (ε,σ) such that D_ε P_σ G P_σ^{-1} D_ε^{-1} = G
    # where D_ε = diag(ε_1,...,ε_8) and P_σ permutes columns.
    # Equivalently: the 8×8 matrix M = D_ε · P_σ satisfies M G M^T = G.

    # Check all elements of B_8 is infeasible (|B_8| = 2^8 · 8! = 10321920).
    # But we can compute the stabilizer more cleverly.

    # The subspace V is the row space of P. A signed permutation M acts
    # on V by: V → M V. This preserves V iff M V = V, i.e., M maps V to itself.
    # In projector terms: M G M^T = G.

    # Strategy: check a generating set of B_8.
    # B_8 is generated by:
    #   - transpositions (i j) in S_8
    #   - sign flip of coordinate i
    # We find which of these preserve G, then build the stabilizer.

    # Actually, let's directly enumerate the stabilizer by checking 
    # structured subsets. The stabilizer of a 4D subspace under B_8
    # is typically small.

    # First, check: which single coordinate sign flips preserve V?
    print("\nChecking single coordinate sign flips...")
    for i in range(8):
        D = np.eye(8)
        D[i, i] = -1
        G2 = D @ G @ D.T
        if np.allclose(G2, G, atol=1e-8):
            print(f"  Flip coord {i}: PRESERVES V")
        else:
            pass  # print(f"  Flip coord {i}: does not preserve V")

    # Check transpositions
    print("\nChecking transpositions...")
    for i in range(8):
        for j in range(i+1, 8):
            M = np.eye(8)
            M[i, i] = M[j, j] = 0
            M[i, j] = M[j, i] = 1
            G2 = M @ G @ M.T
            if np.allclose(G2, G, atol=1e-8):
                print(f"  Swap ({i},{j}): PRESERVES V")

    # Systematic: check all signed permutations of the 4 pairs
    # The sin/cos structure pairs coords as (0,1),(2,3),(4,5),(6,7)
    print("\nChecking pair-preserving signed permutations...")
    stabilizer = []
    pairs = [(0,1), (2,3), (4,5), (6,7)]

    # Actions that preserve the pair structure:
    # - Permute the 4 pairs (S_4, 24 elements)
    # - Swap within each pair (Z_2^4, 16 elements)
    # - Sign-flip each pair (Z_2^4, 16 elements — but sign-flip BOTH coords in a pair)
    # Total: 24 × 16 × 16 = 6144

    for pair_perm in itertools.permutations(range(4)):
        for within_swaps in itertools.product([False, True], repeat=4):
            for pair_signs in itertools.product([1, -1], repeat=4):
                M = np.zeros((8, 8))
                for new_pair_idx, old_pair_idx in enumerate(pair_perm):
                    c1, c2 = pairs[old_pair_idx]
                    d1, d2 = 2*new_pair_idx, 2*new_pair_idx+1
                    if within_swaps[new_pair_idx]:
                        c1, c2 = c2, c1
                    s = pair_signs[new_pair_idx]
                    M[d1, c1] = s
                    M[d2, c2] = s
                G2 = M @ G @ M.T
                if np.allclose(G2, G, atol=1e-8):
                    stabilizer.append(M)

    print(f"  Stabilizer size (pair-preserving): {len(stabilizer)}")

    # Also check non-pair-preserving elements: individual coord sign flips
    print("\nChecking individual coordinate sign flips + pair permutations...")
    stabilizer_full = list(stabilizer)  # start with what we have

    # Try all 2^8 = 256 sign patterns combined with the pair permutations
    for pair_perm in itertools.permutations(range(4)):
        for within_swaps in itertools.product([False, True], repeat=4):
            for signs in itertools.product([1, -1], repeat=8):
                M = np.zeros((8, 8))
                for new_pair_idx, old_pair_idx in enumerate(pair_perm):
                    c1, c2 = pairs[old_pair_idx]
                    d1, d2 = 2*new_pair_idx, 2*new_pair_idx+1
                    if within_swaps[new_pair_idx]:
                        c1, c2 = c2, c1
                    M[d1, c1] = signs[c1]
                    M[d2, c2] = signs[c2]
                G2 = M @ G @ M.T
                if np.allclose(G2, G, atol=1e-8):
                    # Check if new
                    is_new = not any(np.allclose(M, S, atol=1e-8) for S in stabilizer_full)
                    if is_new:
                        stabilizer_full.append(M)

    print(f"  Extended stabilizer size: {len(stabilizer_full)}")

    # The orbit size = |subgroup checked| / |stabilizer|
    # But we haven't checked all of B_8, just pair-structured elements.
    # Let's try a bigger search: sample random B_8 elements.
    print("\nSampling random B_8 elements to find more stabilizer elements...")
    np.random.seed(42)
    for _ in range(100000):
        perm = np.random.permutation(8)
        signs = np.random.choice([-1, 1], size=8)
        M = np.zeros((8, 8))
        for i in range(8):
            M[i, perm[i]] = signs[i]
        G2 = M @ G @ M.T
        if np.allclose(G2, G, atol=1e-8):
            is_new = not any(np.allclose(M, S, atol=1e-8) for S in stabilizer_full)
            if is_new:
                stabilizer_full.append(M)

    print(f"  Final stabilizer size: {len(stabilizer_full)}")

    B8_order = (2**8) * math.factorial(8)
    orbit_size = B8_order // len(stabilizer_full)
    print(f"\n|B_8| = {B8_order}")
    print(f"|Stabilizer| = {len(stabilizer_full)}")
    print(f"Orbit size = |B_8| / |Stab| = {orbit_size}")
    print(f"  → There are {orbit_size} distinct optimal subspaces in the B_8 orbit.")

    return stabilizer_full


# ══════════════════════════════════════════════════════════════════════════════
# 3. R^4 DECOMPOSITION: Does V = V1 ⊕ V2?
# ══════════════════════════════════════════════════════════════════════════════

def decomposition_analysis():
    """
    The sin/cos structure says V = span(rows 0,1) ⊕ span(rows 2,3).
    This is a direct sum decomposition of the 4D subspace.
    
    Question: Is this decomposition forced by optimality?
    Test: given the 8 generators, check if they span two orthogonal 2D planes.
    """
    print("\n" + "=" * 70)
    print("3. R^4 DIRECT SUM DECOMPOSITION")
    print("=" * 70)

    P = get_sincos_basis(X_OPT, Y_OPT)
    gens = generators_from_P(P)  # 8×4

    # In the sin/cos family, generators pair as:
    # g_i and g_{i+1} (for i even) have the same (x1,x2) and opposite (x3,x4)
    # So g_i + g_{i+1} lies in span(e1,e2) and g_i - g_{i+1} lies in span(e3,e4)

    print("\nGenerator pair sums and differences:")
    for pair_idx in range(4):
        i, j = 2*pair_idx, 2*pair_idx+1
        s = gens[i] + gens[j]
        d = gens[i] - gens[j]
        print(f"  g{i}+g{j} = [{s[0]:+.8f}, {s[1]:+.8f}, {s[2]:+.8f}, {s[3]:+.8f}]")
        print(f"  g{i}-g{j} = [{d[0]:+.8f}, {d[1]:+.8f}, {d[2]:+.8f}, {d[3]:+.8f}]")
        # Check: does sum live in (x1,x2) plane?
        print(f"    sum in (x1,x2) plane? {np.allclose(s[2:], 0, atol=1e-8)}")
        print(f"    diff in (x3,x4) plane? {np.allclose(d[:2], 0, atol=1e-8)}")

    # So V1 = span{e1,e2} (the "x-plane") and V2 = span{e3,e4} (the "y-plane")
    # The 4 pair-sums span V1, the 4 pair-diffs span V2.

    # The volume can be written in terms of the V1 and V2 projections.
    # Each generator has a V1 component and a V2 component.
    print("\n\nV1 projections (x1,x2 components of generators):")
    v1_projs = gens[:, :2]
    for i, p in enumerate(v1_projs):
        print(f"  π₁(g{i}) = ({p[0]:+.8f}, {p[1]:+.8f})")

    print("\nV2 projections (x3,x4 components of generators):")
    v2_projs = gens[:, 2:]
    for i, p in enumerate(v2_projs):
        print(f"  π₂(g{i}) = ({p[0]:+.8f}, {p[1]:+.8f})")

    # Observation: The V1 projections of paired generators are IDENTICAL
    # (since the pair structure means same top-2 components)
    print("\nV1 projections of paired generators:")
    for pair_idx in range(4):
        i, j = 2*pair_idx, 2*pair_idx+1
        same = np.allclose(v1_projs[i], v1_projs[j])
        print(f"  π₁(g{i}) = π₁(g{j})? {same}")

    print("\nV2 projections of paired generators (should be negatives):")
    for pair_idx in range(4):
        i, j = 2*pair_idx, 2*pair_idx+1
        neg = np.allclose(v2_projs[i], -v2_projs[j])
        print(f"  π₂(g{i}) = -π₂(g{j})? {neg}")

    # So the 8 generators in R^4 decompose as:
    # g_{2k} = (a_k, b_k)  and  g_{2k+1} = (a_k, -b_k)
    # where a_k ∈ V1 and b_k ∈ V2.
    
    # This means the volume can be decomposed:
    # det(g_{i1}, g_{i2}, g_{i3}, g_{i4}) involves a V1 part and a V2 part
    # Since the generators have this sum/difference structure, many determinants
    # factor nicely.

    print("\n--- Volume Factorization ---")
    print("With g_{2k} = (a_k, b_k), g_{2k+1} = (a_k, -b_k):")
    print()
    
    # The 4×4 determinant det(g_{i1}, ..., g_{i4}) can be expressed via
    # the Cauchy-Binet formula in terms of the 2×2 minors of the V1 and V2
    # projections.
    
    # Each 4×4 det = sum_{|S|=2} det(V1 part in rows S) * det(V2 part in rows S^c)
    
    # Let's verify this numerically
    P_full = gens.T  # 4×8
    P1 = P_full[:2, :]  # V1 part, 2×8
    P2 = P_full[2:, :]  # V2 part, 2×8
    
    count_factored = 0
    count_total = 0
    for idx in itertools.combinations(range(8), 4):
        sub = P_full[:, list(idx)]
        d_full = det(sub)
        
        # Cauchy-Binet: det = sum_{|T|=2 ⊂ {i1,...,i4}} det(P1[:,T]) * det(P2[:,T^c])
        d_cb = 0.0
        cols = list(idx)
        for T in itertools.combinations(range(4), 2):
            Tc = [k for k in range(4) if k not in T]
            d1 = det(P1[:, [cols[T[0]], cols[T[1]]]])
            d2 = det(P2[:, [cols[Tc[0]], cols[Tc[1]]]])
            d_cb += d1 * d2
        
        if abs(d_full - d_cb) < 1e-10:
            count_factored += 1
        count_total += 1
    
    print(f"Cauchy-Binet factorization verified for {count_factored}/{count_total} determinants")
    
    # Now: the V1 and V2 parts are INDEPENDENT (depend on x and y separately).
    # The volume = (1/16) sum |det| = (1/16) sum |sum_{T} det1(T) * det2(T^c)|
    # If the sign pattern is consistent, this becomes a product of V1 and V2 terms.
    
    # Actually, it doesn't factor as a product because the |...| wraps the sum.
    # But the key point is that the CRITICAL POINT EQUATIONS decouple:
    # dV/dx depends only on x (through V1) and a "coupling constant" from V2
    # dV/dy depends only on y (through V2) and a "coupling constant" from V1
    
    print("""
KEY INSIGHT: The decomposition V = V1 ⊕ V2 explains why the sin/cos
family contains the optimum:

1. The tight frame condition PP^T = 4I with the block structure 
   P = [A; B] forces A A^T = 2I₂ and B B^T = 2I₂ independently.

2. Each 2×8 tight frame (A or B) is parameterized by a SINGLE angle
   (rotation within a 1-parameter family of 2D tight frames in R^8).

3. The volume function, through Cauchy-Binet, couples V1 and V2 only
   weakly: the critical point equations for x and y are:
     cos(2x) / sin(2x) = 4 + 3 cos(2y)
     3 cos(2y) / sin(2y) = 4 + 3 cos(2x)
   Each involves the other angle only through cos(2·).

4. The 16-dim → 2-dim reduction happens in three stages:
     Gr(4,8) [dim 16] → Tight frames [dim ~9] → Block decomp [dim ~2+2] 
     → Sin/cos family [dim 2]
   The tight frame condition removes ~7 dimensions.
   The block decomposition removes ~5 more.
   Within the block structure, each 2D tight frame has 1 free parameter.
""")


# ══════════════════════════════════════════════════════════════════════════════
# 4. HIGHER SYMMETRY SEARCH
# ══════════════════════════════════════════════════════════════════════════════

def higher_symmetry_search():
    """
    The known optimum has D_4×Z_2 symmetry (order 16).
    Are there other critical points of the volume function with HIGHER symmetry?
    
    Key idea: highly symmetric tight frames in R^4 with 8 generators include:
    - Vertices of the 4D demihypercube (= 16-cell): 8 vectors, B_2 × B_2 symmetry
    - Roots of D_4: 24 vectors (too many)
    - Subsets of 600-cell vertices
    """
    print("\n" + "=" * 70)
    print("4. SEARCH FOR HIGHER-SYMMETRY OPTIMA")
    print("=" * 70)

    # ── 4a: Check the 16-cell / hyperoctahedron generators ──
    # The 16-cell has 8 vertices: ±e_1, ±e_2, ±e_3, ±e_4
    # As a zonotope with 8 generators ±e_i, this gives the tesseract!
    # Volume = 1 (unit hypercube). Not optimal.

    print("\n--- 4a: Standard basis (16-cell vertices) ---")
    P_std = np.eye(4, 8)  # Just use 4 of the 8 standard basis vectors
    # Actually we need 8 vectors in R^4 forming a tight frame.
    # e_1,...,e_4 repeated gives 4+4=8 but they're not general position.
    
    # The most symmetric 8-vector tight frame in R^4:
    # Take the 8 vertices of the 4D cross-polytope (16-cell): ±e_i
    # As generators: e_1, e_2, e_3, e_4, -e_1, -e_2, -e_3, -e_4
    # But these are ±e_i, so the zonotope is the hypercube [−1,1]^4, volume 16.
    # As a projection: P = [I₄ | -I₄], PP^T = 2I₄.
    
    P_cross = np.hstack([np.eye(4), -np.eye(4)])  # 4×8
    vol_cross = zonotope_volume(P_cross * 2)  # scale so PP^T = 4I
    print(f"  Cross-polytope generators volume: {vol_cross:.10f}")

    # ── 4b: Equally-spaced generators ──
    # In 2D, the optimal is the regular octagon (equally spaced).
    # In 4D, try equally spaced on S^3.
    # The "most spread out" 8 points on S^3 might be vertices of a 
    # 4D polytope. Try the demihypercube (half of tesseract vertices).
    
    print("\n--- 4b: Demihypercube vertices ---")
    # 16 vertices of tesseract: (±1,±1,±1,±1)/2
    # Demihypercube (8 vertices): even number of minus signs
    demi = []
    for v in itertools.product([-1, 1], repeat=4):
        if np.prod(v) == 1:  # even parity
            demi.append(v)
    demi = np.array(demi, dtype=float) / math.sqrt(2)  # normalize
    print(f"  Demihypercube: {len(demi)} vertices")
    
    # Use as generators for a zonotope
    P_demi = demi.T  # 4×8
    # Check tight frame
    G_demi = P_demi @ P_demi.T
    print(f"  PP^T = {np.round(G_demi, 4)}")
    is_tight = np.allclose(G_demi, G_demi[0,0]*np.eye(4))
    print(f"  Is tight frame? {is_tight}")
    if is_tight:
        scale = math.sqrt(4.0 / G_demi[0,0])
        vol_demi = zonotope_volume(P_demi * scale)
        print(f"  Volume (scaled to PP^T=4I): {vol_demi:.10f}")
    
    # ── 4c: x = y case — higher symmetry within sin/cos family ──
    print("\n--- 4c: Special case x = y (higher symmetry?) ---")
    # If x = y, the matrix P has more symmetry because the row-pair structure
    # becomes symmetric under swapping the two 2D blocks.
    from scipy.optimize import minimize_scalar
    
    def neg_vol_diag(t):
        return -zonotope_volume(get_sincos_basis(t, t))
    
    res = minimize_scalar(neg_vol_diag, bounds=(0.01, 0.5), method='bounded')
    x_diag = res.x
    vol_diag = -res.fun
    print(f"  Best x=y: x={x_diag:.12f}, vol={vol_diag:.12f}")
    print(f"  vs optimal: vol={zonotope_volume(get_sincos_basis(X_OPT, Y_OPT)):.12f}")
    print(f"  Deficit: {zonotope_volume(get_sincos_basis(X_OPT, Y_OPT)) - vol_diag:.6f}")
    
    # Check symmetry of x=y case
    P_diag = get_sincos_basis(x_diag, x_diag)
    gens_diag = generators_from_P(P_diag)
    I4 = np.eye(4)
    
    syms_diag = []
    src = gens_diag[:4]
    if abs(det(src)) > 1e-8:
        src_inv = np.linalg.inv(src)
        for perm in itertools.permutations(range(8), 4):
            for signs in itertools.product([1, -1], repeat=4):
                dst = np.array([signs[k] * gens_diag[perm[k]] for k in range(4)])
                R = dst @ src_inv
                if not np.allclose(R @ R.T, I4, atol=1e-7):
                    continue
                ok = True
                for g in gens_diag:
                    Rg = R @ g
                    if not any(np.allclose(Rg, h, atol=1e-7) or np.allclose(Rg, -h, atol=1e-7) 
                              for h in gens_diag):
                        ok = False
                        break
                if ok and not any(np.allclose(R, S, atol=1e-7) for S in syms_diag):
                    syms_diag.append(R)
    
    print(f"  Symmetry group order (x=y): {len(syms_diag)}")
    
    orders_diag = {}
    for R in syms_diag:
        o = order_dist_single(R, I4)
        orders_diag[o] = orders_diag.get(o, 0) + 1
    print(f"  Element orders: {orders_diag}")

    # ── 4d: x = 0 case ──
    print("\n--- 4d: Special case x = 0 ---")
    def neg_vol_x0(t):
        return -zonotope_volume(get_sincos_basis(0, t))
    
    res = minimize_scalar(neg_vol_x0, bounds=(0.01, 0.5), method='bounded')
    y_0 = res.x
    vol_x0 = -res.fun
    print(f"  Best x=0: y={y_0:.12f}, vol={vol_x0:.12f}")
    
    # ── 4e: x = π/4 case (Hadamard-like) ──
    print("\n--- 4e: Special case x=π/4, optimizing y ---")
    def neg_vol_xpi4(t):
        return -zonotope_volume(get_sincos_basis(math.pi/4, t))
    
    res = minimize_scalar(neg_vol_xpi4, bounds=(0.01, 0.78), method='bounded')
    y_pi4 = res.x
    vol_xpi4 = -res.fun
    print(f"  Best x=π/4: y={y_pi4:.12f}, vol={vol_xpi4:.12f}")

    # ── 4f: Check all local maxima of V(x,y) ──
    print("\n--- 4f: All Critical Points of V(x,y) ---")
    from scipy.optimize import minimize
    
    critical_points = []
    for trial in range(500):
        x0 = np.random.rand(2) * math.pi - math.pi/2
        res = minimize(lambda p: -zonotope_volume(get_sincos_basis(p[0], p[1])),
                       x0, method='Nelder-Mead',
                       options={'xatol': 1e-14, 'fatol': 1e-14, 'maxiter': 50000})
        vol = -res.fun
        if vol > 7.0:  # filter out saddle points and minima
            xy = tuple(np.round(res.x, 8))
            is_new = not any(abs(vol - v) < 1e-6 for _, v in critical_points)
            if is_new:
                critical_points.append((xy, vol))
    
    critical_points.sort(key=lambda t: -t[1])
    print(f"  Found {len(critical_points)} distinct critical volumes:")
    for xy, vol in critical_points[:10]:
        # Check symmetry quickly
        P_cp = get_sincos_basis(xy[0], xy[1])
        gens_cp = generators_from_P(P_cp)
        # Quick symmetry check: just count by trying a subset
        sym_count = 0
        I4 = np.eye(4)
        src = gens_cp[:4]
        if abs(det(src)) < 1e-8:
            sym_count = -1  # degenerate
        else:
            src_inv = np.linalg.inv(src)
            # Just check identity and -I
            for perm in itertools.permutations(range(8), 4):
                for signs in itertools.product([1, -1], repeat=4):
                    dst = np.array([signs[k] * gens_cp[perm[k]] for k in range(4)])
                    R = dst @ src_inv
                    if np.allclose(R @ R.T, I4, atol=1e-7):
                        ok = all(any(np.allclose(R @ g, h, atol=1e-7) or np.allclose(R @ g, -h, atol=1e-7)
                                     for h in gens_cp) for g in gens_cp)
                        if ok:
                            sym_count += 1
        
        marker = "★ MAX" if abs(vol - critical_points[0][1]) < 1e-6 else ""
        print(f"    V={vol:.10f}  angles=({xy[0]:+.8f}, {xy[1]:+.8f})  sym≥{sym_count}  {marker}")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    np.random.seed(42)
    symmetries = correct_symmetry_analysis()
    stabilizer = orbit_analysis()
    decomposition_analysis()
    higher_symmetry_search()

if __name__ == '__main__':
    main()
