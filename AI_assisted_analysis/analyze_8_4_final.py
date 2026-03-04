#!/usr/bin/env python3
"""
High-precision analysis: verify the minimal polynomial, compute the volume
symbolically, identify the group, and search for the algebraic structure.
"""

import itertools
import math
import numpy as np
from scipy.linalg import det

def get_bases(x, y):
    sx, cx = math.sin(x), math.cos(x)
    sy, cy = math.sin(y), math.cos(y)
    return np.array([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

def zonotope_volume_direct(x, y):
    P = get_bases(x, y)
    vol = 0.0
    for idx in itertools.combinations(range(8), 4):
        vol += abs(det(P[:, list(idx)]))
    return vol / 16.0

# ── Part 1: Verify the symbolic volume formula ──────────────────────────────

def verify_symbolic_formula():
    """Check that V(x,y) = 2 + sin(2x)/2 + 3sin(2y)/2 + 2cos(2x) + 2cos(2y) + 3cos(2x)cos(2y)/2
    matches the explicit computation at several points."""
    
    print("=" * 70)
    print("VERIFYING SYMBOLIC VOLUME FORMULA")
    print("=" * 70)
    
    def V_formula(x, y):
        return (2 + math.sin(2*x)/2 + 3*math.sin(2*y)/2 
                + 2*math.cos(2*x) + 2*math.cos(2*y) 
                + 3*math.cos(2*x)*math.cos(2*y)/2)
    
    test_points = [
        (0.07348, 0.20329),  # near optimum
        (0.05, 0.15),
        (0.01, 0.30),
        (0.1, 0.1),
        (math.pi/16, math.pi/8),
    ]
    
    for x, y in test_points:
        v_explicit = zonotope_volume_direct(x, y)
        v_formula = V_formula(x, y)
        diff = abs(v_explicit - v_formula)
        status = "✓" if diff < 1e-8 else "✗"
        print(f"  ({x:.5f}, {y:.5f}): explicit={v_explicit:.12f}, formula={v_formula:.12f}, diff={diff:.2e} {status}")


# ── Part 2: High-precision optimization ──────────────────────────────────────

def high_precision_optimization():
    """Use mpmath for high-precision optimization."""
    print("\n" + "=" * 70)
    print("HIGH-PRECISION ANALYSIS")
    print("=" * 70)
    
    try:
        import mpmath
        mpmath.mp.dps = 50  # 50 decimal places
        
        def V_mp(x, y):
            s2x = mpmath.sin(2*x)
            c2x = mpmath.cos(2*x)
            s2y = mpmath.sin(2*y)
            c2y = mpmath.cos(2*y)
            return 2 + s2x/2 + 3*s2y/2 + 2*c2x + 2*c2y + 3*c2x*c2y/2
        
        def dVdx_mp(x, y):
            return mpmath.cos(2*x) - 4*mpmath.sin(2*x) - 3*mpmath.sin(2*x)*mpmath.cos(2*y)
        
        def dVdy_mp(x, y):
            return 3*mpmath.cos(2*y) - 4*mpmath.sin(2*y) - 3*mpmath.cos(2*x)*mpmath.sin(2*y)
        
        # Solve the system dV/dx = 0, dV/dy = 0
        from mpmath import findroot
        
        x0, y0 = mpmath.mpf('0.07348'), mpmath.mpf('0.2033')
        
        sol = findroot(lambda x, y: (dVdx_mp(x, y), dVdy_mp(x, y)), (x0, y0))
        x_opt, y_opt = sol
        
        print(f"Optimal x = {x_opt}")
        print(f"Optimal y = {y_opt}")
        
        u = mpmath.cos(2*x_opt)
        v = mpmath.cos(2*y_opt)
        V = V_mp(x_opt, y_opt)
        
        print(f"\ncos(2x) = {u}")
        print(f"cos(2y) = {v}")
        print(f"Volume  = {V}")
        
        # Verify the minimal polynomial
        poly_val = (2025*u**8 + 10800*u**7 + 24840*u**6 + 19848*u**5 
                    - 18368*u**4 - 51744*u**3 - 24448*u**2 + 21504*u + 16384)
        print(f"\nMinimal polynomial at u: {poly_val}")
        
        # Try PSLQ for volume
        print("\n--- PSLQ for Volume ---")
        # Try to express V as a linear combination of algebraic numbers
        targets = [mpmath.mpf(1), mpmath.sqrt(2), mpmath.sqrt(3), mpmath.sqrt(5), 
                   mpmath.sqrt(6), mpmath.sqrt(7), mpmath.sqrt(10), mpmath.sqrt(15)]
        
        # V = c0 + c1*sqrt(2) + c2*sqrt(3) + ...
        vec = [V] + targets
        try:
            rel = mpmath.pslq(vec)
            if rel:
                print(f"  PSLQ relation found: {rel}")
                # Reconstruct: rel[0]*V + rel[1]*1 + rel[2]*sqrt(2) + ... = 0
                # V = -(rel[1] + rel[2]*sqrt(2) + ...)/rel[0]
        except:
            print("  No PSLQ relation found with sqrt(2,3,5,6,7,10,15)")
        
        # Try to find minimal polynomial of V using PSLQ
        print("\n--- Minimal polynomial of V ---")
        powers = [V**k for k in range(13)]
        rel = mpmath.pslq(powers)
        if rel:
            print(f"  Minimal polynomial coefficients (V^0 to V^{len(rel)-1}): {rel}")
            # Verify
            check = sum(c * V**k for k, c in enumerate(rel))
            print(f"  Verification: sum = {check}")
        else:
            print("  No integer relation found up to degree 12")
        
        # Also try 16*V
        print("\n--- Minimal polynomial of 16*V ---")
        W = 16*V
        powers_W = [W**k for k in range(13)]
        rel_W = mpmath.pslq(powers_W)
        if rel_W:
            print(f"  Minimal polynomial coefficients: {rel_W}")
            check = sum(c * W**k for k, c in enumerate(rel_W))
            print(f"  Verification: sum = {check}")
        
        # Try V^2
        print("\n--- Minimal polynomial of V^2 ---")
        V2 = V**2
        powers_V2 = [V2**k for k in range(9)]
        rel_V2 = mpmath.pslq(powers_V2)
        if rel_V2:
            print(f"  Minimal polynomial coefficients: {rel_V2}")
            check = sum(c * V2**k for k, c in enumerate(rel_V2))
            print(f"  Verification: sum = {check}")
        
        return x_opt, y_opt, V, u, v
        
    except ImportError:
        print("mpmath not available")
        return None


# ── Part 3: Verify and factor the minimal polynomial ────────────────────────

def analyze_polynomial():
    """Study the degree 8 polynomial for cos(2x)."""
    print("\n" + "=" * 70)
    print("MINIMAL POLYNOMIAL ANALYSIS")
    print("=" * 70)
    
    import sympy
    from sympy import Poly, Symbol, factor, roots, Rational, sqrt, nroots
    
    u = Symbol('u')
    p = 2025*u**8 + 10800*u**7 + 24840*u**6 + 19848*u**5 - 18368*u**4 - 51744*u**3 - 24448*u**2 + 21504*u + 16384
    
    poly = Poly(p, u)
    print(f"Polynomial: {p}")
    print(f"Degree: {poly.degree()}")
    
    # Check irreducibility
    factors = sympy.factor_list(p)
    print(f"Factor list: {factors}")
    
    if len(factors[1]) == 1 and factors[1][0][1] == 1:
        print("The polynomial is IRREDUCIBLE over Q.")
    
    # Find all roots numerically
    numerical_roots = nroots(p)
    print(f"\nNumerical roots:")
    for i, r in enumerate(numerical_roots):
        # Check if root corresponds to a valid cos(2x) value
        r_float = complex(r)
        is_real = abs(r_float.imag) < 1e-10
        in_range = is_real and -1 <= r_float.real <= 1
        status = "✓ valid cos value" if in_range else ""
        print(f"  root {i}: {r}  {status}")
    
    # Coefficients
    print(f"\nCoefficients: {poly.all_coeffs()}")
    print(f"Leading coeff: 2025 = 45^2 = (3^2 * 5)^2")
    print(f"Constant term: 16384 = 2^14")
    print(f"Discriminant (large number): ...")
    
    # Galois group analysis (may be slow)
    # print("\nGalois group...")
    # try:
    #     from sympy.polys.numberfields import galois_group
    #     G = galois_group(poly)
    #     print(f"Galois group: {G}")
    # except:
    #     print("Could not compute Galois group")
    
    # Check the derivative of V with respect to u at the root
    # From V = 2 + u/(2(4+3v)) + 9v/(2(4+3u)) + 2u + 2v + 3uv/2
    # where v satisfies its own polynomial...
    
    # Let me also find the polynomial for v = cos(2y)
    v = Symbol('v')
    eq1 = u**2 * (1 + (4 + 3*v)**2) - (4 + 3*v)**2
    eq2 = v**2 * (9 + (4 + 3*u)**2) - (4 + 3*u)**2
    
    # Resultant to eliminate u
    from sympy import resultant
    res_v = resultant(Poly(sympy.expand(eq1), u), Poly(sympy.expand(eq2), u), u)
    res_v_factored = sympy.factor(res_v)
    print(f"\nResultant for v: {res_v_factored}")
    
    # The polynomial for v should be related to the one for u by the symmetry of the system
    

# ── Part 4: Group identification ─────────────────────────────────────────────

def identify_group():
    """Identify the symmetry group of order 16."""
    print("\n" + "=" * 70)
    print("SYMMETRY GROUP IDENTIFICATION")
    print("=" * 70)
    
    # The 16 symmetries from earlier analysis.
    # Generators: a = S8, b = S15
    # a = 90° rotation in (x1,x2) and negate x3
    # b = swap(x1,x2) and swap(x3,x4)
    # Relations: a^4 = b^2 = 1, ab = -ba (where - is central inversion)
    
    S8 = np.array([[0,-1,0,0],[1,0,0,0],[0,0,-1,0],[0,0,0,1]], dtype=float)
    S15 = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]], dtype=float)
    I4 = np.eye(4)
    
    a = S8
    b = S15
    
    # Verify relations
    a2 = a @ a
    a4 = a2 @ a2
    b2 = b @ b
    
    print(f"a^2 = {a2.tolist()}")
    print(f"a^4 = {a4.tolist()} (= I? {np.allclose(a4, I4)})")
    print(f"b^2 = {b2.tolist()} (= I? {np.allclose(b2, I4)})")
    
    ab = a @ b
    ba = b @ a
    neg_ba = -ba
    print(f"\nab = {ab.tolist()}")
    print(f"-ba = {neg_ba.tolist()}")
    print(f"ab = -ba? {np.allclose(ab, neg_ba)}")
    
    # Generate all group elements
    elements = [I4]
    queue = [a, b]
    
    while queue:
        g = queue.pop(0)
        is_new = True
        for e in elements:
            if np.allclose(g, e, atol=1e-10):
                is_new = False
                break
        if is_new:
            elements.append(g)
            for h in elements[:-1]:
                queue.append(g @ h)
                queue.append(h @ g)
    
    print(f"\nGroup order: {len(elements)}")
    
    # Compute element orders
    order_counts = {}
    for g in elements:
        power = I4.copy()
        for k in range(1, 17):
            power = power @ g
            if np.allclose(power, I4, atol=1e-10):
                order_counts[k] = order_counts.get(k, 0) + 1
                break
    
    print(f"Element order distribution: {order_counts}")
    
    # Compute center
    center = []
    for g in elements:
        is_central = True
        for h in elements:
            if not np.allclose(g @ h, h @ g, atol=1e-10):
                is_central = False
                break
        if is_central:
            center.append(g)
    
    print(f"Center order: {len(center)}")
    
    # The group with order 16, center order 4, element orders {1:1, 2:7, 4:8}
    # or similar, can help identify it.
    
    # Known groups of order 16 and their properties:
    # D_8: orders {1:1, 2:5, 4:2, 8:8}, center Z_2 - NO (center too small)
    # Q_16: orders {1:1, 4:6, 8:8, 16:1}, center Z_2 - NO
    # D_4 × Z_2: orders {1:1, 2:7, 4:8}, center Z_2 × Z_2 - POSSIBLE!
    # Q_8 × Z_2: orders {1:1, 2:3, 4:12}, center Z_2 × Z_2 - POSSIBLE
    # Z_4 × Z_4: orders {1:1, 2:3, 4:12}, center = G - NO (not abelian)
    # (Z_4 × Z_2) ⋊ Z_2: varies
    
    print(f"\nComparing with known groups of order 16:")
    print(f"  D_4 × Z_2: orders {{1:1, 2:7, 4:8}}, center Z_2×Z_2")
    print(f"  Q_8 × Z_2: orders {{1:1, 2:3, 4:12}}, center Z_2×Z_2")
    print(f"  Our group:  orders {order_counts}, center order {len(center)}")
    
    if order_counts.get(1,0) == 1 and order_counts.get(2,0) == 7 and order_counts.get(4,0) == 8:
        print(f"\n  ★ The symmetry group is D_4 × Z_2 (dihedral group of square × Z_2)")
    elif order_counts.get(1,0) == 1 and order_counts.get(2,0) == 3 and order_counts.get(4,0) == 12:
        print(f"\n  ★ The symmetry group is Q_8 × Z_2 (quaternion group × Z_2)")


# ── Part 5: Literature connection ────────────────────────────────────────────

def literature_connection():
    """Discuss connections to known mathematical problems."""
    print("\n" + "=" * 70)
    print("MATHEMATICAL CONTEXT AND SIGNIFICANCE")
    print("=" * 70)
    
    print("""
1. PROBLEM DEFINITION:
   Find the 4D subspace W ⊂ R^8 that maximizes Vol_4(π_W([0,1]^8))
   where π_W is the orthogonal projection onto W.

   Equivalently: maximize the zonotope volume generated by the projections
   of the 8 standard basis vectors of R^8 onto a 4D subspace.

2. KNOWN RESULTS FOR 2n → n:
   ┌─────────┬──────────────────────┬───────────────┬──────────────────┐
   │ n       │ Shadow polytope      │ Volume        │ Symmetry order   │
   ├─────────┼──────────────────────┼───────────────┼──────────────────┤
   │ 1 (2→1) │ Line segment         │ √2            │ 2                │
   │ 2 (4→2) │ Regular octagon      │ 2+2√2 ≈ 4.83 │ 16 (D_8)         │
   │ 3 (6→3) │ Rhombic triaconta-   │ 5√5 ≈ 11.18  │ 120 (I_h)        │
   │         │ hedron               │               │                  │
   │ 4 (8→4) │ Equilateral zonotope │ ≈ 7.845       │ 16               │
   │         │ Z(8,4)               │               │                  │
   └─────────┴──────────────────────┴───────────────┴──────────────────┘

3. BREAK IN THE PATTERN:
   For n=2: related to the silver ratio (1+√2), D_8 symmetry
   For n=3: related to the golden ratio φ = (1+√5)/2, icosahedral symmetry
   For n=4: NO connection to any known special constant or root system.
   
   The optimal angles satisfy transcendental equations:
     cot(2x) = 4 + 3cos(2y)
     3cot(2y) = 4 + 3cos(2x)
   
   cos(2x) satisfies an IRREDUCIBLE degree-8 polynomial over Q:
     2025u^8 + 10800u^7 + 24840u^6 + 19848u^5 
     - 18368u^4 - 51744u^3 - 24448u^2 + 21504u + 16384 = 0
   
   This means the optimal angles are NOT rational multiples of π and
   cannot be expressed using radicals of small degree.

4. THE 4D POLYTOPE:
   • f-vector: (128, 352, 336, 112) — maximum for Z(8,4)
   • All 352 edges have equal length (equilateral zonotope)
   • 112 facets are parallelotopes (combinatorial cubes) of 7 types
   • 336 ridges (2-faces) are all parallelograms
   • Symmetry group of order 16
   • NOT a regular polytope, NOT related to any root system
   
5. WHY THE PATTERN BREAKS:
   The 4→2 and 6→3 cases benefit from the existence of highly symmetric
   arrangements of vectors (the silver ratio for D_4, the icosahedral 
   arrangement for I_h). These are related to the fact that:
   - In 2D, the regular n-gon is the unique (up to rotation) arrangement
     of n unit vectors maximizing the zonotope area.
   - In 3D, the icosahedral arrangement of 6 vectors is optimal.
   
   In 4D, there is no analogous highly symmetric arrangement of 8 vectors 
   that maximizes the zonotope volume. The optimal arrangement has only 
   modest symmetry (order 16 out of a maximum possible 384 for the 
   hyperoctahedral group B_4).

6. CONNECTIONS TO OTHER PROBLEMS:
   The shadow volume problem is related to:
   - Geometric tomography (Bourgain-Lindenstrauss-Milman)
   - Randomized dimension reduction (Johnson-Lindenstrauss lemma)
   - The Busemann-Petty problem for zonotopes
   - Coding theory (finding good projections for error correction)
   
   The specific value ≈ 7.845 and the degree-8 minimal polynomial suggest 
   this is a "generic" optimization problem without special algebraic 
   structure, unlike the n=2 and n=3 cases.
""")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    verify_symbolic_formula()
    result = high_precision_optimization()
    analyze_polynomial()
    identify_group()
    literature_connection()

if __name__ == '__main__':
    main()
