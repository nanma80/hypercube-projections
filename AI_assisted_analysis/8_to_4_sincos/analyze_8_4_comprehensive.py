#!/usr/bin/env python3
"""
Comprehensive analysis of the 8D -> 4D maximum-shadow projection.

The projection matrix P(x,y) is parameterized by two angles x, y:
  Row 0: [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx]
  Row 1: [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx]
  Row 2: [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy]
  Row 3: [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]

where sx=sin(x), cx=cos(x), sy=sin(y), cy=cos(y).
Rows are orthogonal with norm 2, so P P^T = 4 I_4.
Volume = (1/16) sum_{|S|=4} |det(P_S)|.
"""

import itertools
import math
import numpy as np
from scipy.linalg import orth, det
from scipy.optimize import minimize_scalar, minimize
from scipy.spatial import ConvexHull
from collections import Counter

# ── Part 0: Projection matrix and volume ─────────────────────────────────────

def get_bases(x, y):
    sx, cx = math.sin(x), math.cos(x)
    sy, cy = math.sin(y), math.cos(y)
    return np.array([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

def zonotope_volume(bases):
    """Volume = sum |det(Q_S)| where Q is orthonormal basis for subspace."""
    n = bases.shape[1]
    k = bases.shape[0]
    Q = orth(bases.T)
    vol = 0.0
    for idx in itertools.combinations(range(n), k):
        vol += abs(det(Q[list(idx), :]))
    return vol

def zonotope_volume_direct(x, y):
    """Volume using the known structure P P^T = 4I, so Q = P^T/2."""
    P = get_bases(x, y)
    n, k = P.shape[1], P.shape[0]
    vol = 0.0
    for idx in itertools.combinations(range(n), k):
        sub = P[:, list(idx)]
        vol += abs(det(sub))
    return vol / 16.0

def get_cube_vertices(dim):
    return np.array(list(itertools.product([-0.5, 0.5], repeat=dim)))


# ── Part 1: Numerical optimization ──────────────────────────────────────────

def optimize():
    """Find optimal angles by grid search + local optimization."""
    best_vol = 0
    best_xy = (0, 0)
    # Grid search
    for ix in range(1, 50):
        for iy in range(1, 50):
            x = ix * 0.002
            y = iy * 0.01
            v = zonotope_volume_direct(x, y)
            if v > best_vol:
                best_vol = v
                best_xy = (x, y)
    # Local optimization
    res = minimize(lambda p: -zonotope_volume_direct(p[0], p[1]),
                   best_xy, method='Nelder-Mead',
                   options={'xatol': 1e-14, 'fatol': 1e-14, 'maxiter': 100000})
    x_opt, y_opt = res.x
    vol_opt = -res.fun
    return x_opt, y_opt, vol_opt


# ── Part 2: Symbolic analysis of determinants ────────────────────────────────

def symbolic_determinants():
    """Compute all 70 determinants symbolically and group by expression."""
    import sympy
    from sympy import symbols, sin, cos, Matrix, trigsimp, factor, expand_trig

    x, y = symbols('x y', real=True, positive=True)
    sx, cx = sin(x), cos(x)
    sy, cy = sin(y), cos(y)

    P = Matrix([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

    results = []
    for idx in itertools.combinations(range(8), 4):
        sub = P[:, list(idx)]
        d = sub.det()
        d_simple = trigsimp(d)
        results.append((idx, d_simple))

    # Group by simplified expression (ignoring sign)
    from collections import defaultdict
    groups = defaultdict(list)
    for idx, d in results:
        key = trigsimp(d**2)
        groups[str(key)].append((idx, d))

    return x, y, results, groups


def symbolic_volume_and_gradient():
    """Build the volume function (with signs fixed at the optimum) and differentiate."""
    import sympy
    from sympy import symbols, sin, cos, Matrix, trigsimp, diff, Rational, solve, Eq

    x, y = symbols('x y', real=True, positive=True)
    sx, cx = sin(x), cos(x)
    sy, cy = sin(y), cos(y)

    P = Matrix([
        [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
        [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
        [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
        [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

    # Compute all determinants
    dets = []
    for idx in itertools.combinations(range(8), 4):
        sub = P[:, list(idx)]
        d = sub.det()
        dets.append((idx, trigsimp(d)))

    # Evaluate signs at the numerical optimum
    x_num, y_num = 0.07348072421223484, 0.20328550404073642
    signed_volume = sympy.Integer(0)
    for idx, d in dets:
        val = float(d.subs([(x, x_num), (y, y_num)]))
        sign = 1 if val >= 0 else -1
        signed_volume += sign * d

    signed_volume = trigsimp(signed_volume / 16)
    print("\n=== Symbolic volume (with signs from optimum) ===")
    print(f"V(x,y) = {signed_volume}")

    # Gradient
    dVdx = trigsimp(diff(signed_volume, x))
    dVdy = trigsimp(diff(signed_volume, y))
    print(f"\ndV/dx = {dVdx}")
    print(f"\ndV/dy = {dVdy}")

    return signed_volume, dVdx, dVdy, x, y


# ── Part 3: Convex hull and f-vector ─────────────────────────────────────────

def analyze_hull(x, y):
    """Compute the convex hull and extract the f-vector."""
    P = get_bases(x, y)
    Q = P.T / 2  # 8×4 orthonormal-ish (actually P^T/2 since PP^T=4I)
    vertices = get_cube_vertices(8)  # 256 × 8
    projected = vertices @ Q  # 256 × 4

    hull = ConvexHull(projected)
    hull_vert_coords = projected[hull.vertices]

    # Unique facet equations (3-cells)
    eps = 1e-10
    rounded_eq = np.round(hull.equations, 10)
    unique_eq = np.unique(rounded_eq, axis=0)
    n_cells = len(unique_eq)

    # For each cell, find its vertices
    cells = []
    for eq in unique_eq:
        on_facet = []
        for vi in range(len(projected)):
            val = np.dot(projected[vi], eq[:4]) + eq[4]
            if val > -eps:
                on_facet.append(vi)
        cells.append(tuple(sorted(on_facet)))

    # Faces (2-faces): each cell is a 3D polytope; find its 2-faces
    # by intersecting pairs of cells
    cell_vertex_sets = [set(c) for c in cells]
    faces = set()
    edges = set()
    all_hull_verts = set(hull.vertices)

    # Find 2-faces as intersections of pairs of 3-cells sharing >= 3 vertices
    for i in range(len(cells)):
        for j in range(i+1, len(cells)):
            shared = cell_vertex_sets[i] & cell_vertex_sets[j]
            if len(shared) >= 3:
                face = tuple(sorted(shared))
                faces.add(face)

    # Find edges as intersections of 2-faces sharing >= 2 vertices
    face_list = list(faces)
    face_vertex_sets = [set(f) for f in face_list]
    for i in range(len(face_list)):
        for j in range(i+1, len(face_list)):
            shared = face_vertex_sets[i] & face_vertex_sets[j]
            if len(shared) == 2:
                edge = tuple(sorted(shared))
                edges.add(edge)

    n_vertices = len(hull.vertices)
    n_edges = len(edges)
    n_faces = len(faces)

    print(f"\n=== f-vector of the 4D polytope ===")
    print(f"Vertices: {n_vertices}")
    print(f"Edges:    {n_edges}")
    print(f"2-faces:  {n_faces}")
    print(f"3-cells:  {n_cells}")
    print(f"Euler characteristic (V - E + F - C): {n_vertices - n_edges + n_faces - n_cells}")
    print(f"(Should be 0 for a 4D polytope)")

    # Analyze cell types (by vertex count)
    cell_sizes = Counter(len(c) for c in cells)
    print(f"\n3-cell types by vertex count: {dict(cell_sizes)}")

    # Analyze face types
    face_sizes = Counter(len(f) for f in faces)
    print(f"2-face types by vertex count: {dict(face_sizes)}")

    # Vertex distances
    dists = np.linalg.norm(hull_vert_coords, axis=1)
    unique_dists = np.unique(np.round(dists, 8))
    print(f"\nDistinct vertex distances from origin: {unique_dists}")
    dist_counts = Counter(np.round(dists, 8))
    print(f"Vertex distance multiplicities: {dict(sorted(dist_counts.items()))}")

    # Edge lengths
    edge_lengths = []
    for e in edges:
        v0, v1 = list(e)
        length = np.linalg.norm(projected[v0] - projected[v1])
        edge_lengths.append(length)
    unique_edge_lengths = np.unique(np.round(edge_lengths, 8))
    print(f"\nDistinct edge lengths: {unique_edge_lengths}")
    el_counts = Counter(np.round(edge_lengths, 8))
    print(f"Edge length multiplicities: {dict(sorted(el_counts.items()))}")

    return hull, projected, cells, faces, edges, hull_vert_coords


# ── Part 4: Symmetry group ──────────────────────────────────────────────────

def find_symmetry_group(hull_vert_coords, eps=1e-7):
    """Find all orthogonal transformations mapping vertex set to itself.
    
    Strategy: enumerate signed permutation matrices (hyperoctahedral group B_4,
    order 384) and check each one. If needed, also try other transformations.
    """
    n = len(hull_vert_coords)
    vert_set = set()
    for v in hull_vert_coords:
        vert_set.add(tuple(np.round(v, 7)))

    symmetries = []

    # Check all signed permutations (B_4, order 384)
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product([-1, 1], repeat=4):
            T = np.zeros((4, 4))
            for i in range(4):
                T[i, perm[i]] = signs[i]
            # Apply T to all vertices
            transformed = hull_vert_coords @ T.T
            trans_set = set()
            for v in transformed:
                trans_set.add(tuple(np.round(v, 7)))
            if trans_set == vert_set:
                symmetries.append(T)

    print(f"\n=== Symmetry analysis ===")
    print(f"Signed permutation symmetries found: {len(symmetries)}")
    print(f"(B_4 = hyperoctahedral group has order {math.factorial(4) * 2**4})")

    # Check if we need to look beyond signed permutations
    # Try some rotation matrices that aren't signed permutations
    extra_count = 0
    # Try rotations by common angles in each plane
    angles = [math.pi/6, math.pi/4, math.pi/3, math.pi/5, math.pi/8, math.pi/10, math.pi/12]
    for theta in angles:
        for phi in angles:
            # Rotation in (x1,x2) plane by theta and (x3,x4) plane by phi
            ct, st = math.cos(theta), math.sin(theta)
            cp, sp = math.cos(phi), math.sin(phi)
            T = np.array([
                [ct, -st, 0, 0],
                [st,  ct, 0, 0],
                [0, 0, cp, -sp],
                [0, 0, sp,  cp]
            ])
            transformed = hull_vert_coords @ T.T
            trans_set = set(tuple(np.round(v, 7)) for v in transformed)
            if trans_set == vert_set:
                # Check if already found
                already = False
                for S in symmetries:
                    if np.allclose(S, T, atol=eps):
                        already = True
                        break
                if not already:
                    extra_count += 1
                    symmetries.append(T)
                    print(f"  Extra symmetry found: rot({theta*180/math.pi:.1f}°, {phi*180/math.pi:.1f}°)")

    if extra_count == 0:
        print("No additional rotation symmetries found beyond signed permutations.")

    # Print the symmetry matrices
    if len(symmetries) <= 20:
        print("\nSymmetry matrices:")
        for i, T in enumerate(symmetries):
            print(f"  S{i}: {T.tolist()}")

    return symmetries


# ── Part 5: Try to identify the volume as a closed form ──────────────────────

def identify_volume(vol):
    """Try to match the volume against closed-form expressions."""
    print(f"\n=== Volume identification ===")
    print(f"Numerical volume: {vol:.15f}")
    print(f"16 * vol = {16*vol:.15f}")

    phi = (1 + math.sqrt(5)) / 2
    sqrt2 = math.sqrt(2)
    sqrt3 = math.sqrt(3)
    sqrt5 = math.sqrt(5)

    candidates = {
        "8*sqrt(5-1)/2": 8 * (sqrt5 - 1) / 2,
        "4*phi^2": 4 * phi**2,
        "8/phi + 4*phi": 8/phi + 4*phi,
        "2*(1+sqrt(2))^2": 2 * (1+sqrt2)**2,
        "4*(1+sqrt(3))": 4 * (1+sqrt3),
        "4 + 4*sqrt(2) + sqrt(5)": 4 + 4*sqrt2 + sqrt5,
        "8*sqrt(2)*phi/sqrt(3)": 8*sqrt2*phi/sqrt3,
        "16*phi/sqrt(5)": 16*phi/sqrt5,
        "3 + 5*sqrt(2)": 3 + 5*sqrt2,
        "1 + 5*sqrt(2) + sqrt(5)": 1 + 5*sqrt2 + sqrt5,
        "2^(5/2) + pi": 2**2.5 + math.pi,
        "5*sqrt(5)/2 + sqrt(2)": 5*sqrt5/2 + sqrt2,
        "7 + sqrt(2)/sqrt(3)": 7 + sqrt2/sqrt3,
        "sqrt(61.5)": math.sqrt(61.5),
        "4*sqrt(2) + sqrt(3)": 4*sqrt2 + sqrt3,
        "3*sqrt(2) + 2*sqrt(5)": 3*sqrt2 + 2*sqrt5,
        "2*sqrt(2) + 3*sqrt(3)": 2*sqrt2 + 3*sqrt3,
    }

    print("\nCandidate closed forms:")
    for name, val in sorted(candidates.items(), key=lambda kv: abs(kv[1] - vol)):
        diff = val - vol
        print(f"  {name:40s} = {val:.15f}  diff = {diff:+.2e}")

    # Also try: is vol^2 nice? is 16*vol nice? etc.
    print(f"\nvol^2 = {vol**2:.15f}")
    print(f"vol^3 = {vol**3:.15f}")
    print(f"(16*vol)^2 = {(16*vol)**2:.15f}")

    # Try minimal polynomial approach using sympy
    try:
        import sympy
        from sympy import nsimplify, Rational
        approx = sympy.nsimplify(vol, rational=False, tolerance=1e-8)
        print(f"\nsympy nsimplify(vol): {approx} = {float(approx):.15f}")
        approx16 = sympy.nsimplify(16*vol, rational=False, tolerance=1e-8)
        print(f"sympy nsimplify(16*vol): {approx16} = {float(approx16):.15f}")
    except Exception as e:
        print(f"nsimplify failed: {e}")


# ── Part 6: Identify the optimal angles ─────────────────────────────────────

def identify_angles(x, y):
    """Try to identify what the optimal angles are."""
    print(f"\n=== Angle identification ===")
    print(f"x = {x:.15f} rad = {math.degrees(x):.12f}°")
    print(f"y = {y:.15f} rad = {math.degrees(y):.12f}°")
    print(f"x + y = {x+y:.15f} rad = {math.degrees(x+y):.12f}°")
    print(f"x - y = {x-y:.15f} rad = {math.degrees(x-y):.12f}°")
    print(f"2x = {2*x:.15f} rad = {math.degrees(2*x):.12f}°")
    print(f"2y = {2*y:.15f} rad = {math.degrees(2*y):.12f}°")
    print(f"2x + 2y = {2*(x+y):.15f} rad = {math.degrees(2*(x+y)):.12f}°")

    print(f"\nTrig values at optimal angles:")
    print(f"sin(x) = {math.sin(x):.15f}")
    print(f"cos(x) = {math.cos(x):.15f}")
    print(f"sin(y) = {math.sin(y):.15f}")
    print(f"cos(y) = {math.cos(y):.15f}")
    print(f"sin(2x) = {math.sin(2*x):.15f}")
    print(f"cos(2x) = {math.cos(2*x):.15f}")
    print(f"sin(2y) = {math.sin(2*y):.15f}")
    print(f"cos(2y) = {math.cos(2*y):.15f}")
    print(f"tan(x) = {math.tan(x):.15f}")
    print(f"tan(y) = {math.tan(y):.15f}")
    print(f"tan(x)*tan(y) = {math.tan(x)*math.tan(y):.15f}")
    print(f"tan(x)/tan(y) = {math.tan(x)/math.tan(y):.15f}")
    print(f"sin(x)/sin(y) = {math.sin(x)/math.sin(y):.15f}")
    print(f"cos(x)/cos(y) = {math.cos(x)/math.cos(y):.15f}")

    # Try nsimplify on trig values
    import sympy
    for name, val in [("sin(x)", math.sin(x)), ("cos(x)", math.cos(x)),
                      ("sin(y)", math.sin(y)), ("cos(y)", math.cos(y)),
                      ("tan(x)", math.tan(x)), ("tan(y)", math.tan(y)),
                      ("sin(2x)", math.sin(2*x)), ("cos(2x)", math.cos(2*x)),
                      ("sin(2y)", math.sin(2*y)), ("cos(2y)", math.cos(2*y)),
                      ("x+y", x+y), ("x-y", x-y), ("x/y", x/y),
                      ("(x+y)/pi", (x+y)/math.pi), ("(x-y)/pi", (x-y)/math.pi),
                      ("x/pi", x/math.pi), ("y/pi", y/math.pi),
                      ]:
        try:
            approx = sympy.nsimplify(val, rational=False, tolerance=1e-8)
            print(f"  nsimplify({name}) = {approx}")
        except:
            pass


# ── Part 7: Analyze the generator zonotope structure ─────────────────────────

def analyze_generators(x, y):
    """Analyze the 8 generators of the zonotope."""
    P = get_bases(x, y)
    gens = P.T / 2  # 8 generators, each a 4D vector

    print(f"\n=== Generator analysis ===")
    print("8 generators (columns of P/2):")
    for i in range(8):
        g = gens[i]
        print(f"  g{i+1} = [{', '.join(f'{v:+.8f}' for v in g)}]  |g| = {np.linalg.norm(g):.10f}")

    # Pairwise inner products
    print("\nPairwise inner products (g_i · g_j):")
    for i, j in itertools.combinations(range(8), 2):
        ip = np.dot(gens[i], gens[j])
        if abs(ip) > 1e-10:
            print(f"  g{i+1}·g{j+1} = {ip:+.10f}")

    # Pairwise angles
    print("\nPairwise angles between generators:")
    angle_vals = []
    for i, j in itertools.combinations(range(8), 2):
        cos_angle = np.dot(gens[i], gens[j]) / (np.linalg.norm(gens[i]) * np.linalg.norm(gens[j]))
        cos_angle = np.clip(cos_angle, -1, 1)
        angle = math.degrees(math.acos(abs(cos_angle)))
        angle_vals.append((i+1, j+1, angle, cos_angle))
    angle_vals.sort(key=lambda t: t[2])
    for i, j, angle, cos_a in angle_vals:
        print(f"  g{i}-g{j}: {angle:8.4f}° (cos = {cos_a:+.8f})")


# ── Part 8: Detailed 3-cell analysis ────────────────────────────────────────

def analyze_cells_detail(cells, projected, eps=1e-8):
    """Analyze the 3-cells in detail: are they all the same type?"""
    print(f"\n=== Detailed 3-cell analysis ===")

    cell_types = {}
    for ci, cell in enumerate(cells):
        verts = projected[list(cell)]
        center = np.mean(verts, axis=0)
        rel = verts - center

        # Compute pairwise distances
        dists = []
        for i in range(len(rel)):
            for j in range(i+1, len(rel)):
                dists.append(np.linalg.norm(rel[i] - rel[j]))
        dists = sorted(np.round(dists, 8))
        key = tuple(dists)

        if key not in cell_types:
            cell_types[key] = []
        cell_types[key].append(ci)

    print(f"Number of distinct cell types: {len(cell_types)}")
    for k, (key, indices) in enumerate(cell_types.items()):
        print(f"\n  Type {k+1}: {len(indices)} cells, {len(cells[indices[0]])} vertices each")
        print(f"    Pairwise distances: {list(key)}")
        # Check if this is a known polytope
        n_dists = len(key)
        unique_dists = sorted(set(key))
        dist_counts = Counter(key)
        print(f"    Unique distances: {unique_dists}")
        print(f"    Distance multiplicities: {dict(sorted(dist_counts.items()))}")

        # Check central symmetry
        cell = cells[indices[0]]
        verts = projected[list(cell)]
        center = np.mean(verts, axis=0)
        rel = verts - center
        is_centrally_sym = True
        for v in rel:
            if not any(np.linalg.norm(v + w) < eps for w in rel):
                is_centrally_sym = False
                break
        print(f"    Centrally symmetric: {is_centrally_sym}")

    return cell_types


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("COMPREHENSIVE ANALYSIS: 8D -> 4D Max-Shadow Projection")
    print("=" * 70)

    # Step 1: Optimize
    print("\n--- Step 1: Numerical optimization ---")
    x_opt, y_opt, vol_opt = optimize()
    print(f"Optimal angles: x = {x_opt:.15f}, y = {y_opt:.15f}")
    print(f"Optimal volume: {vol_opt:.15f}")
    print(f"Target volume:  7.84468782041")

    # Verify with known values from v3
    x_known, y_known = 0.07348072421223484, 0.20328550404073642
    vol_known = zonotope_volume_direct(x_known, y_known)
    print(f"\nUsing known angles: vol = {vol_known:.15f}")

    # Use the more precise values
    x, y = x_known, y_known
    vol = vol_known

    # Step 2: Analyze generators
    analyze_generators(x, y)

    # Step 3: Identify angles
    identify_angles(x, y)

    # Step 4: Hull and f-vector
    print("\n--- Step 4: Convex hull analysis ---")
    hull, projected, cells, faces, edges, hull_verts = analyze_hull(x, y)

    # Step 5: Detailed cell analysis
    analyze_cells_detail(cells, projected)

    # Step 6: Symmetry
    print("\n--- Step 6: Symmetry group ---")
    symmetries = find_symmetry_group(hull_verts)

    # Step 7: Volume identification
    identify_volume(vol)

    # Step 8: Symbolic analysis (may be slow)
    print("\n--- Step 8: Symbolic analysis ---")
    try:
        V_sym, dVdx, dVdy, xs, ys = symbolic_volume_and_gradient()
    except Exception as e:
        print(f"Symbolic analysis failed: {e}")


if __name__ == '__main__':
    main()
