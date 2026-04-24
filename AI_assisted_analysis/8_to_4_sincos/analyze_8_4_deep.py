#!/usr/bin/env python3
"""
Deeper analysis: derive critical point equations, find the symmetry group
more carefully, and try to identify the polytope.
"""

import itertools
import math
import numpy as np
from scipy.linalg import det
from scipy.spatial import ConvexHull
from scipy.optimize import fsolve

# ── Projection matrix ────────────────────────────────────────────────────────

def get_bases(x, y):
    sx, cx = math.sin(x), math.cos(x)
    sy, cy = math.sin(y), math.cos(y)
    return np.array([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

def get_cube_vertices(dim):
    return np.array(list(itertools.product([-0.5, 0.5], repeat=dim)))

# ── Part 1: Critical point equations ─────────────────────────────────────────

def critical_point_analysis():
    """
    From the symbolic volume:
    V(x,y) = sin(2x)/2 + 3sin(2y)/2 + 2cos(2x) + 2cos(2y) + 3cos(2x)cos(2y)/2 + 2
    
    Setting gradient to zero gives:
    (I)  cot(2x) = 4 + 3cos(2y)
    (II) 3cot(2y) = 4 + 3cos(2x)
    
    Let u = cos(2x), v = cos(2y), then sin(2x) = sqrt(1-u^2), sin(2y) = sqrt(1-v^2).
    (I)  u/sqrt(1-u^2) = 4 + 3v
    (II) 3v/sqrt(1-v^2) = 4 + 3u
    """
    print("=" * 70)
    print("CRITICAL POINT EQUATIONS")
    print("=" * 70)
    
    # Verify
    x, y = 0.07348072421223484, 0.20328550404073642
    u, v = math.cos(2*x), math.cos(2*y)
    s2x, s2y = math.sin(2*x), math.sin(2*y)
    
    print(f"\nNumerical values:")
    print(f"  x = {x:.15f}, y = {y:.15f}")
    print(f"  u = cos(2x) = {u:.15f}")
    print(f"  v = cos(2y) = {v:.15f}")
    print(f"  sin(2x) = {s2x:.15f}")
    print(f"  sin(2y) = {s2y:.15f}")
    
    lhs1 = u / s2x
    rhs1 = 4 + 3*v
    print(f"\n  Eq(I): cot(2x) = {lhs1:.15f}, 4+3cos(2y) = {rhs1:.15f}, diff = {lhs1-rhs1:.2e}")
    
    lhs2 = 3*v / s2y
    rhs2 = 4 + 3*u
    print(f"  Eq(II): 3cot(2y) = {lhs2:.15f}, 4+3cos(2x) = {rhs2:.15f}, diff = {lhs2-rhs2:.2e}")
    
    # The system in u,v:
    # u^2 = (4+3v)^2 (1-u^2)  =>  u^2 = (4+3v)^2 / (1 + (4+3v)^2)
    # 9v^2 = (4+3u)^2 (1-v^2) =>  v^2 = (4+3u)^2 / (9 + (4+3u)^2)
    
    # Substitute the first into the second.
    # From first: u = (4+3v)/sqrt(1+(4+3v)^2)  (positive root)
    # Then: v^2 = (4+3u)^2 / (9 + (4+3u)^2)
    
    # This defines v as a function of itself, h(v) = v.
    
    def equations(v_val):
        u_val = (4 + 3*v_val) / math.sqrt(1 + (4 + 3*v_val)**2)
        v_pred = (4 + 3*u_val) / math.sqrt(9 + (4 + 3*u_val)**2)
        return v_pred - v_val
    
    # Solve for v
    from scipy.optimize import brentq
    v_sol = brentq(equations, 0.5, 0.999)
    u_sol = (4 + 3*v_sol) / math.sqrt(1 + (4 + 3*v_sol)**2)
    
    print(f"\n  Solved: u = {u_sol:.15f}, v = {v_sol:.15f}")
    print(f"  Check: u_known = {u:.15f}, v_known = {v:.15f}")
    
    # Now let's find the minimal polynomial of u (and v)
    # Eliminate v from the system:
    # u^2(1+(4+3v)^2) = (4+3v)^2   =>  (4+3v)^2 = u^2/(1-u^2)
    # =>  4+3v = u/sqrt(1-u^2)  =>  v = (u/sqrt(1-u^2) - 4)/3
    #
    # Substitute into second equation:
    # 9v^2(1-v^2)^(-1) = (4+3u)^2 ... wait, let me use:
    # 9v^2 = (4+3u)^2(1-v^2)
    # 9v^2 + (4+3u)^2 v^2 = (4+3u)^2
    # v^2(9 + (4+3u)^2) = (4+3u)^2
    # v = (4+3u)/sqrt(9+(4+3u)^2)
    
    # And from above: v = (u/sqrt(1-u^2) - 4)/3
    # So: (u/sqrt(1-u^2) - 4)/3 = (4+3u)/sqrt(9+(4+3u)^2)
    
    # Let w = sqrt(1-u^2) (= sin(2x)). Then u^2 + w^2 = 1.
    # v = (u/w - 4)/3 = (u - 4w)/(3w)
    # Also: (4+3u)/sqrt(9+(4+3u)^2)
    
    # (u - 4w)/(3w) = (4+3u)/sqrt(9+(4+3u)^2)
    # (u - 4w) * sqrt(9+(4+3u)^2) = 3w(4+3u)
    
    # Square both sides:
    # (u-4w)^2 (9+(4+3u)^2) = 9w^2(4+3u)^2
    
    # Expand (u-4w)^2 = u^2 - 8uw + 16w^2 = u^2 - 8uw + 16(1-u^2) = 16 - 15u^2 - 8uw
    # Actually let me keep w^2 = 1-u^2.
    
    # (u-4w)^2(9+(4+3u)^2) = 9(1-u^2)(4+3u)^2
    
    # Let me expand everything in terms of u and w, then substitute w^2 = 1-u^2.
    
    # This is getting very messy. Let me use sympy.
    print("\n--- Using sympy to find minimal polynomial ---")
    import sympy
    from sympy import symbols, sqrt, solve, Eq, Rational, simplify, expand, factor, Poly, resultant
    
    u_sym, v_sym = symbols('u v', positive=True)
    
    # System: u^2(1+(4+3v)^2) = (4+3v)^2
    #         v^2(9+(4+3u)^2) = (4+3u)^2
    
    eq1 = u_sym**2 * (1 + (4 + 3*v_sym)**2) - (4 + 3*v_sym)**2
    eq2 = v_sym**2 * (9 + (4 + 3*u_sym)**2) - (4 + 3*u_sym)**2
    
    print(f"  eq1 = {eq1}")
    print(f"  eq2 = {eq2}")
    
    # Compute resultant to eliminate v
    print("  Computing resultant (this may take a moment)...")
    eq1_poly = Poly(expand(eq1), v_sym)
    eq2_poly = Poly(expand(eq2), v_sym)
    
    print(f"  eq1 degree in v: {eq1_poly.degree()}")
    print(f"  eq2 degree in v: {eq2_poly.degree()}")
    
    res = resultant(eq1_poly, eq2_poly, v_sym)
    res_simplified = factor(res)
    
    print(f"\n  Resultant (factored): {res_simplified}")
    
    # Find the minimal polynomial of u
    res_poly = Poly(res_simplified, u_sym)
    print(f"  Degree of resultant in u: {res_poly.degree()}")
    
    # Evaluate at u = u_sol to verify
    res_factors = sympy.factor_list(res)
    print(f"\n  Factor list: {res_factors}")
    
    # Check which factor vanishes at u = u_sol
    for i, (fac, mult) in enumerate(res_factors[1]):
        val = float(fac.subs(u_sym, u_sol))
        print(f"  Factor {i} (mult {mult}): {fac} = {val:.2e} at u={u_sol:.10f}")
    
    return u_sol, v_sol


# ── Part 2: Comprehensive symmetry check ─────────────────────────────────────

def comprehensive_symmetry_check():
    """
    Check all possible symmetries more carefully.
    Use all 256 projected vertices (not just hull vertices) for robustness.
    """
    print("\n" + "=" * 70)
    print("COMPREHENSIVE SYMMETRY ANALYSIS")
    print("=" * 70)
    
    x, y = 0.07348072421223484, 0.20328550404073642
    P = get_bases(x, y)
    Q = P.T / 2  # 8×4
    vertices = get_cube_vertices(8)
    projected = vertices @ Q  # 256×4
    
    hull = ConvexHull(projected)
    hull_vert_indices = set(hull.vertices)
    hull_verts = projected[hull.vertices]
    n_hull = len(hull.vertices)
    
    print(f"Total projected vertices: {len(projected)}")
    print(f"Hull vertices: {n_hull}")
    
    # Build a set of all projected points (rounded)
    decimals = 9
    all_set = set(tuple(np.round(v, decimals)) for v in projected)
    hull_set = set(tuple(np.round(v, decimals)) for v in hull_verts)
    
    # Part A: Check all 384 signed permutations
    sym_count = 0
    symmetries = []
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product([-1, 1], repeat=4):
            T = np.zeros((4, 4))
            for i in range(4):
                T[i, perm[i]] = signs[i]
            transformed = hull_verts @ T.T
            trans_set = set(tuple(np.round(v, decimals)) for v in transformed)
            if trans_set == hull_set:
                sym_count += 1
                symmetries.append(T)
    
    print(f"\nSigned permutation symmetries: {sym_count}")
    
    # Part B: Check generator-based symmetries
    # The generators are g_1,...,g_8 (rows of Q = P^T/2)
    gens = Q  # 8×4
    
    # For each orthogonal T that permutes generators (up to sign),
    # we can enumerate: for each possible mapping of g_1 to ±g_j,
    # check if the resulting T extends to a valid symmetry.
    
    # More practical: check T that we derived analytically
    analytical_transforms = [
        # (x1,x2,x3,x4) → (x2, -x1, x3, -x4): 90° rot in (12) plane, negate x4
        np.array([[0,1,0,0],[-1,0,0,0],[0,0,1,0],[0,0,0,-1]], dtype=float),
        # (x1,x2,x3,x4) → (-x2, x1, -x3, x4): opposite 90° rot in (12), negate x3
        np.array([[0,-1,0,0],[1,0,0,0],[0,0,-1,0],[0,0,0,1]], dtype=float),
        # (x1,x2,x3,x4) → (-x2, x1, x3, -x4)
        np.array([[0,-1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,-1]], dtype=float),
        # (x1,x2,x3,x4) → (x2, -x1, -x3, x4)
        np.array([[0,1,0,0],[-1,0,0,0],[0,0,-1,0],[0,0,0,1]], dtype=float),
    ]
    
    print("\nChecking analytically-derived transforms:")
    for T in analytical_transforms:
        # Check how T acts on generators
        gen_map = []
        valid = True
        for i in range(8):
            tg = T @ gens[i]
            found = False
            for j in range(8):
                if np.allclose(tg, gens[j], atol=1e-10):
                    gen_map.append((i+1, j+1, '+'))
                    found = True
                    break
                elif np.allclose(tg, -gens[j], atol=1e-10):
                    gen_map.append((i+1, j+1, '-'))
                    found = True
                    break
            if not found:
                valid = False
                break
        
        if valid:
            # Also verify on hull vertices
            transformed = hull_verts @ T.T
            trans_set = set(tuple(np.round(v, decimals)) for v in transformed)
            is_sym = trans_set == hull_set
            print(f"  T = {T.tolist()}")
            print(f"    Generator mapping: {gen_map}")
            print(f"    Is symmetry: {is_sym}")
            if is_sym and not any(np.allclose(T, S) for S in symmetries):
                symmetries.append(T)
                sym_count += 1
        else:
            tg = T @ gens[0]
            print(f"  T = {T.tolist()} -> NOT a valid generator permutation")
            print(f"    T(g1) = {tg}, no match found")
    
    # Part C: Search for symmetries by trying to map g1 to each possible ±g_j
    # and building the 4×4 matrix from 4 independent constraints.
    print("\nSearching for symmetries by building T from generator maps...")
    
    extra_syms = []
    # Try all possible images of g1, g3 (they span different directions)
    for j1 in range(8):
        for s1 in [1, -1]:
            for j3 in range(8):
                for s3 in [1, -1]:
                    if j1 == j3:
                        continue
                    # T maps g1 → s1*g_j1, g3 → s3*g_j3
                    # g1 and g3 are linearly independent in R^4
                    # Need two more constraints. Use g5 and g7.
                    for j5 in range(8):
                        for s5 in [1, -1]:
                            for j7 in range(8):
                                for s7 in [1, -1]:
                                    if len({j1,j3,j5,j7}) < 4:
                                        # The target generators should be independent
                                        pass  # might still work
                                    
                                    # Build system: T @ [g1|g3|g5|g7] = [s1*g_j1|s3*g_j3|s5*g_j5|s7*g_j7]
                                    A = np.column_stack([gens[0], gens[2], gens[4], gens[6]])
                                    B = np.column_stack([s1*gens[j1], s3*gens[j3], s5*gens[j5], s7*gens[j7]])
                                    
                                    d = np.linalg.det(A)
                                    if abs(d) < 1e-10:
                                        continue
                                    
                                    T = B @ np.linalg.inv(A)
                                    
                                    # Check orthogonality
                                    if not np.allclose(T @ T.T, np.eye(4), atol=1e-8):
                                        continue
                                    
                                    # Check all 8 generators
                                    all_ok = True
                                    for i in range(8):
                                        tg = T @ gens[i]
                                        if not any(np.allclose(tg, s*gens[k], atol=1e-8)
                                                   for k in range(8) for s in [1,-1]):
                                            all_ok = False
                                            break
                                    
                                    if all_ok:
                                        # Check it's not already found
                                        already = any(np.allclose(T, S, atol=1e-8) for S in symmetries)
                                        if not already:
                                            # Verify on hull
                                            transformed = hull_verts @ T.T
                                            trans_set = set(tuple(np.round(v, decimals)) for v in transformed)
                                            if trans_set == hull_set:
                                                symmetries.append(T)
                                                extra_syms.append(T)
    
    print(f"Additional symmetries found: {len(extra_syms)}")
    print(f"Total symmetry group order: {len(symmetries)}")
    
    for T in extra_syms:
        print(f"  New symmetry: {np.round(T, 8).tolist()}")
    
    # Verify group closure
    print("\nVerifying group closure...")
    products = set()
    for S1 in symmetries:
        for S2 in symmetries:
            prod = S1 @ S2
            found = False
            for S in symmetries:
                if np.allclose(prod, S, atol=1e-8):
                    found = True
                    break
            if not found:
                key = tuple(np.round(prod.flatten(), 7))
                products.add(key)
    
    if products:
        print(f"  Group NOT closed! Missing {len(products)} products.")
        # Add missing products and iterate
        for key in products:
            T = np.array(key).reshape(4,4)
            symmetries.append(T)
        print(f"  After adding: {len(symmetries)} elements")
    else:
        print(f"  Group is closed. Order = {len(symmetries)}")
    
    return symmetries


# ── Part 3: Volume at the critical point ─────────────────────────────────────

def volume_at_critical_point():
    """
    V = 2 + sin(2x)/2 + 3sin(2y)/2 + 2cos(2x) + 2cos(2y) + 3cos(2x)cos(2y)/2
    
    Using cot(2x) = 4 + 3cos(2y) and 3cot(2y) = 4 + 3cos(2x):
    sin(2x) = cos(2x)/(4+3cos(2y)) = u/(4+3v)
    sin(2y) = 3cos(2y)/(4+3cos(2x)) = 3v/(4+3u)
    """
    print("\n" + "=" * 70)
    print("VOLUME AT CRITICAL POINT")
    print("=" * 70)
    
    x, y = 0.07348072421223484, 0.20328550404073642
    u = math.cos(2*x)
    v = math.cos(2*y)
    
    V = 2 + math.sin(2*x)/2 + 3*math.sin(2*y)/2 + 2*u + 2*v + 3*u*v/2
    print(f"  V = {V:.15f}")
    
    # Substitute sin(2x) = u/(4+3v), sin(2y) = 3v/(4+3u):
    V2 = 2 + u/(2*(4+3*v)) + 9*v/(2*(4+3*u)) + 2*u + 2*v + 3*u*v/2
    print(f"  V (substituted) = {V2:.15f}")
    
    # Express in terms of u only (using v = v(u))
    # From: v = (4+3u)/sqrt(9+(4+3u)^2)  [derived from eq II]
    B = 4 + 3*u
    v_from_u = B / math.sqrt(9 + B**2)
    print(f"  v from u: {v_from_u:.15f} vs {v:.15f}")
    
    # And: u = (4+3v)/sqrt(1+(4+3v)^2)  [from eq I]
    A = 4 + 3*v
    u_from_v = A / math.sqrt(1 + A**2)
    print(f"  u from v: {u_from_v:.15f} vs {u:.15f}")
    
    # Let's try to find V as a function of u alone
    v_u = (4+3*u) / math.sqrt(9 + (4+3*u)**2)
    sin2x = u / (4 + 3*v_u)
    sin2y = 3*v_u / (4 + 3*u)
    
    V3 = 2 + sin2x/2 + 3*sin2y/2 + 2*u + 2*v_u + 3*u*v_u/2
    print(f"  V (from u only) = {V3:.15f}")
    
    # Simplify sin2y when v = (4+3u)/sqrt(9+(4+3u)^2):
    # sin2y = 3v/(4+3u) = 3(4+3u)/(sqrt(9+(4+3u)^2)*(4+3u)) = 3/sqrt(9+(4+3u)^2)
    sin2y_simple = 3 / math.sqrt(9 + (4+3*u)**2)
    print(f"  sin(2y) simplified = {sin2y_simple:.15f} vs {math.sin(2*y):.15f}")
    
    # And 4+3v = 4 + 3(4+3u)/sqrt(9+(4+3u)^2) = (4*sqrt(9+B^2) + 3B)/sqrt(9+B^2)
    # where B = 4+3u
    # sin2x = u/(4+3v) = u*sqrt(9+B^2)/(4*sqrt(9+B^2)+3B)
    
    # This is getting complicated. Let me try with sympy.
    print("\n--- Using sympy to simplify V(u) ---")
    import sympy
    from sympy import symbols, sqrt as Sqrt, Rational, simplify, nsimplify
    
    u_s = symbols('u', positive=True)
    B_s = 4 + 3*u_s
    v_s = B_s / Sqrt(9 + B_s**2)
    A_s = 4 + 3*v_s
    sin2x_s = u_s / A_s
    sin2y_s = 3 / Sqrt(9 + B_s**2)  # simplified form
    
    V_s = 2 + sin2x_s/2 + Rational(3,2)*sin2y_s + 2*u_s + 2*v_s + Rational(3,2)*u_s*v_s
    V_simplified = simplify(V_s)
    print(f"  V(u) = {V_simplified}")
    
    # Also compute numerically for some test u values
    print(f"\n  V({float(u):.10f}) = {float(V_simplified.subs(u_s, u)):.15f}")
    
    # Try nsimplify on the volume
    vol_approx = nsimplify(V, rational=False, tolerance=1e-6)
    print(f"\n  nsimplify(V, tol=1e-6) = {vol_approx} ≈ {float(vol_approx):.15f}")
    vol_approx2 = nsimplify(V, [sympy.sqrt(2), sympy.sqrt(3), sympy.sqrt(5), sympy.pi], 
                             rational=False, tolerance=1e-4)
    print(f"  nsimplify(V, constants, tol=1e-4) = {vol_approx2} ≈ {float(vol_approx2):.15f}")
    
    # Try: is V^2 nice?
    V2_approx = nsimplify(V**2, rational=False, tolerance=1e-4)
    print(f"  nsimplify(V^2, tol=1e-4) = {V2_approx} ≈ {float(V2_approx):.15f}, actual V^2 = {V**2:.15f}")
    
    # Try 16*V
    V16_approx = nsimplify(16*V, rational=False, tolerance=1e-4)
    print(f"  nsimplify(16*V, tol=1e-4) = {V16_approx} ≈ {float(V16_approx):.15f}")


# ── Part 4: Zonotope as known polytope ───────────────────────────────────────

def identify_polytope():
    """Compare f-vector and structure with known 4-polytopes."""
    print("\n" + "=" * 70)
    print("POLYTOPE IDENTIFICATION")
    print("=" * 70)
    
    print("""
The polytope has:
  f-vector: (128, 352, 336, 112)
  All edges equal length (1/√2)
  All 3-cells are combinatorial cubes (8 vertices, 12 edges each)
  All 2-faces are quadrilaterals (parallelograms)
  7 distinct types of 3-cells (by metric)
  Each group has 16 cells
  Symmetry group order: ≥16

Known zonotopes from 2n → n projections:
  4→2: Regular octagon (f-vector: 8, 8; symmetry D_8, order 16)
  6→3: Rhombic triacontahedron (f-vector: 32, 60, 30; symmetry I_h, order 120)
  8→4: Our polytope (f-vector: 128, 352, 336, 112; symmetry order ??)

For comparison, the 8→4 zonotope with maximum symmetry would use generators
related to the D_4 or F_4 root system. The 24-cell has f-vector (24,96,96,24).

Our polytope has 112 = 7×16 cells, 336 = 7×48 faces, 352 edges, 128 = 2^7 vertices.
Note: 128 = 2^7, and a zonotope from 8 generators in R^4 has at most
  2 * C(8-1, 4-1) = 2 * C(7,3) = 2*35 = 70 faces in general position,
  but with special structure can have more.

Actually, for a zonotope generated by n vectors in R^d, the number of vertices
is at most 2*sum_{k=0}^{d-1} C(n-1,k). For n=8, d=4:
  2*(C(7,0)+C(7,1)+C(7,2)+C(7,3)) = 2*(1+7+21+35) = 2*64 = 128.
  
This matches! So 128 is the MAXIMUM number of vertices for a zonotope from 
8 generators in R^4. This means the generators are in "general position"
(no 4 generators are coplanar).

Number of 2-faces of a zonotope = sum over pairs of generators that are 
NOT coplanar of 1... actually for a zonotope:
  # faces (zones) = C(n,2) for general position = C(8,2) = 28 zones, but
  each zone consists of multiple faces.
  
For a zonotope: #cells = 2*C(n-1, d-1) = 2*C(7,3) = 70... hmm that's not 112.
Wait, that formula is for simple zonotopes (all vertices have minimal valence).
""")
    
    # Count combinatorial data
    x, y = 0.07348072421223484, 0.20328550404073642
    P = get_bases(x, y)
    Q = P.T / 2
    gens = Q  # 8×4, rows are generators
    
    # Check if any 4 generators are coplanar (linearly dependent)
    print("Checking if generators are in general position...")
    for idx in itertools.combinations(range(8), 4):
        sub = gens[list(idx)]
        d = abs(det(sub))
        if d < 1e-10:
            print(f"  DEGENERATE: generators {idx} have det = {d:.2e}")
    
    # Count zones: for each pair (i,j), the 2D plane spanned by g_i, g_j
    # determines a "zone" of faces
    print(f"\nNumber of generator pairs: {len(list(itertools.combinations(range(8), 2)))}")
    
    # The zonotope formula for n generators in R^d, general position:
    # f_k = 2 * C(n, d-k) * ... this is the Zaslavsky formula
    # Actually the f-vector of a zonotope Z(n,d) with n generators in general position in R^d:
    # f_k = 2 * sum_{j=k}^{d} (-1)^{j-k} C(j,k) * C(n,j) ... no.
    # 
    # For a zonotope with n generators in R^d in general position:
    # The number of k-faces is 2 * C(n-1, d-k-1) * C(n, k)  / C(d, k) ... 
    # This doesn't look right either.
    
    # Let me just look up the formula.
    # f_0 (vertices) = 2 * sum_{i=0}^{d-1} C(n-1, i)
    # For n=8, d=4: 2*(1+7+21+35) = 128 ✓
    
    # f_{d-1} (facets) = 2 * C(n-1, d-1) for simple zonotopes... but 2*C(7,3) = 70 ≠ 112
    
    # So our zonotope is NOT simple. There are degenerate vertices where more
    # than d facets meet.
    
    # Actually, the formula for the number of facets of the zonotope Z(V) where
    # V = {v_1,...,v_n} in R^d is:
    # Each facet corresponds to choosing (d-1) generators that span a hyperplane,
    # and the two sides. So f_{d-1} = 2 * (# hyperplane arrangements from choosing d-1 generators)
    
    # For general position: f_{d-1} = 2 * C(n, d-1) = 2 * C(8,3) = 2*56 = 112. ✓!
    
    print(f"\nGeneral position formula for facets: 2*C(8,3) = {2*56}")
    print(f"Actual facets: 112 ✓")
    
    # f_1 (edges): For general position zonotope, each edge is parallel to one generator.
    # For generator g_i, the edges parallel to it form a "zone".
    # The number of edges parallel to g_i = number of vertices of the projection
    # of the zonotope onto the hyperplane perpendicular to g_i, which is a 
    # (d-1)-dimensional zonotope with (n-1) generators.
    # 
    # For d=4, n=8: the number of edges parallel to one generator is the number of
    # vertices of a 3D zonotope with 7 generators = 2*(C(6,0)+C(6,1)+C(6,2)) = 2*(1+6+15) = 44.
    # Total edges = 8 * 44 / 2 ... no, each edge is counted once.
    # Actually, each edge is parallel to exactly one generator, so
    # total edges = sum over generators of (# edges parallel to that generator)
    # = 8 * 44 = 352. ✓!
    # Wait, 8 * 44 = 352. Each edge is counted once since each edge is parallel 
    # to exactly one generator. Yes!
    
    print(f"General position formula for edges: 8 * 2*(1+6+15) = 8*44 = {8*44}")
    print(f"Actual edges: 352 ✓")
    
    # f_2 (2-faces): Each 2-face is a parallelogram, determined by a pair of generators.
    # For each pair (g_i, g_j), the number of 2-faces in the corresponding zone is
    # related to the number of vertices of the zonotope projected onto the (d-2)-dim
    # space perpendicular to span(g_i, g_j).
    # For d=4, n=8: the projection is a 2D zonotope with 6 generators, which has
    # 2*C(5,0) + 2*C(5,1) = 2+10 = 12 vertices.
    # Total 2-faces = C(8,2) * 12 = 28 * 12 = 336. ✓!
    
    print(f"General position formula for 2-faces: C(8,2)*2*(1+5) = 28*12 = {28*12}")
    print(f"Actual 2-faces: 336 ✓")
    
    print(f"""
CONCLUSION: The f-vector (128, 352, 336, 112) matches perfectly with a 
GENERAL POSITION zonotope from 8 generators in R^4.

The polytope is a 4-dimensional zonotope Z(8,4), which is the Minkowski sum 
of 8 line segments in 4D. Equivalently, it is the shadow (projection) of 
the 8-dimensional unit cube onto a 4-dimensional subspace.

All edges are equal (length 1/√2 ≈ 0.7071), confirming it is an 
EQUILATERAL zonotope (isozonohedron generalized to 4D).

The 7 types of 3-cells (each appearing 16 times, total 112) are all 
combinatorial cubes (parallelotopes) — consistent with zonotope theory 
where all facets are zonotopes of one lower dimension.
""")

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    u_sol, v_sol = critical_point_analysis()
    comprehensive_symmetry_check()
    volume_at_critical_point()
    identify_polytope()

if __name__ == '__main__':
    main()
