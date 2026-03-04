#!/usr/bin/env python3
"""
Generate plots and document generators for the 10→3 and 8→3 max-volume
projections of hypercubes.

10→3: Rhombic enneacontahedron (I_h symmetry, order 120)
      Generators = 10 three-fold axes of the icosahedron
      Zero free parameters — uniquely determined by icosahedral symmetry.

8→3:  D_7 × Z_2 zonotope (order 28)
      Generators = 7 vectors on a cone + 1 axial vector
      One free parameter θ (cone half-angle).
"""

import itertools
import math
import numpy as np
from scipy.linalg import det
from scipy.spatial import ConvexHull
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

PHI = (1 + math.sqrt(5)) / 2  # golden ratio


# ══════════════════════════════════════════════════════════════════════════════
# 10→3: Rhombic Enneacontahedron
# ══════════════════════════════════════════════════════════════════════════════

def build_10_3_generators():
    """
    The 10 generators are the three-fold rotation axes of the icosahedron,
    i.e., the directions through the face centers of the 20 triangular faces
    (10 antipodal pairs → 10 directions).

    Icosahedron vertices: cyclic permutations of (0, ±1, ±φ).
    Face center = (v_a + v_b + v_c) / 3 for each triangular face.
    Normalized to unit vectors, they automatically form a tight frame:
      PP^T = (10/3) I_3.

    SYMBOLIC FORM (exact):
      Each generator is (v_a + v_b + v_c) / |v_a + v_b + v_c|
      where v_a, v_b, v_c are vertices of an icosahedral face.
      All entries involve only 1, φ = (1+√5)/2, and their sums.
    """
    ico = []
    for cyc in [(0,1,2), (1,2,0), (2,0,1)]:
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = [0.0, 0.0, 0.0]
                v[cyc[1]] = s1
                v[cyc[2]] = s2 * PHI
                ico.append(v)
    ico = np.array(ico)

    # Find triangular faces (edge length = 2)
    faces = []
    for i in range(12):
        for j in range(i+1, 12):
            for k in range(j+1, 12):
                ds = [np.linalg.norm(ico[a]-ico[b]) for a,b in [(i,j),(j,k),(i,k)]]
                if all(abs(d - 2.0) < 0.01 for d in ds):
                    faces.append((i,j,k))

    # Face centers as unit vectors
    fc = np.array([(ico[f[0]]+ico[f[1]]+ico[f[2]])/3 for f in faces])
    fc_unit = fc / np.linalg.norm(fc, axis=1, keepdims=True)

    # Pick one from each antipodal pair
    axes = []
    used = set()
    for i in range(len(fc_unit)):
        if i in used:
            continue
        for j in range(i+1, len(fc_unit)):
            if j not in used and np.allclose(fc_unit[i], -fc_unit[j], atol=1e-6):
                axes.append(fc_unit[i])
                used.add(i)
                used.add(j)
                break
    return np.array(axes), ico, faces


# ══════════════════════════════════════════════════════════════════════════════
# 8→3: D_7 × Z_2 zonotope
# ══════════════════════════════════════════════════════════════════════════════

def build_8_3_generators(theta=1.065776577865147):
    """
    The 8 generators split into 7 + 1:
      g_k = a·(sinθ·cos(2πk/7), sinθ·sin(2πk/7), cosθ)   for k = 0,...,6
      g_7 = (0, 0, b)

    The tight-frame condition PP^T = (8/3)I determines a and b from θ:
      a² = 16 / (21·sin²θ)
      b² = 8·(1 - 3·cos²θ) / (3·sin²θ)

    Valid for cos²θ < 1/3 (i.e., θ > arccos(1/√3) ≈ 0.9553).

    Optimal: θ ≈ 1.06578, a ≈ 0.9974, b ≈ 1.0182.
    """
    st, ct = math.sin(theta), math.cos(theta)
    a = math.sqrt(16 / (21 * st**2))
    b = math.sqrt(8 * (1 - 3*ct**2) / (3 * st**2))

    gens = []
    for k in range(7):
        phi = 2 * math.pi * k / 7
        gens.append([a*st*math.cos(phi), a*st*math.sin(phi), a*ct])
    gens.append([0, 0, b])
    return np.array(gens), a, b


def build_zonotope(gens, n):
    """Build zonotope vertices as projections of [0,1]^n cube vertices."""
    cube = np.array(list(itertools.product([0, 1], repeat=n)))
    proj = cube @ gens
    proj -= proj.mean(axis=0)
    return proj


def get_zonotope_faces(proj):
    """Extract parallelogram faces from a zonotope convex hull.
    
    Merges coplanar triangle pairs from ConvexHull into quadrilaterals.
    """
    hull = ConvexHull(proj)
    
    # Group triangles by their outward normal (coplanar faces share a normal)
    normals = hull.equations[:, :3]  # outward normals
    offsets = hull.equations[:, 3]
    
    # Merge coplanar triangles: two triangles sharing an edge and the same normal
    used = [False] * len(hull.simplices)
    quads = []
    
    for i in range(len(hull.simplices)):
        if used[i]:
            continue
        for j in range(i+1, len(hull.simplices)):
            if used[j]:
                continue
            # Check if coplanar (same normal and offset)
            if (np.allclose(normals[i], normals[j], atol=1e-6) and
                abs(offsets[i] - offsets[j]) < 1e-6):
                # Find shared edge and unique vertices
                si = set(hull.simplices[i])
                sj = set(hull.simplices[j])
                shared = si & sj
                if len(shared) == 2:
                    unique_i = (si - shared).pop()
                    unique_j = (sj - shared).pop()
                    shared = list(shared)
                    # Order as a proper quadrilateral: unique_i, shared[0], unique_j, shared[1]
                    quad = [unique_i, shared[0], unique_j, shared[1]]
                    quads.append(quad)
                    used[i] = True
                    used[j] = True
                    break
    
    # Any remaining unpaired triangles (shouldn't happen for zonotopes)
    tris = [hull.simplices[i] for i in range(len(hull.simplices)) if not used[i]]
    
    return quads, tris


def plot_zonotope(ax, proj, title, facecolor, edgecolor):
    """Plot a zonotope with proper quadrilateral faces."""
    quads, tris = get_zonotope_faces(proj)
    
    # Draw quadrilateral faces
    for quad in quads:
        verts = proj[quad]
        poly = Poly3DCollection([verts], alpha=1.0, facecolor=facecolor,
                                edgecolor=edgecolor, linewidth=0.5)
        ax.add_collection3d(poly)
    
    # Draw any remaining triangles
    for tri in tris:
        verts = proj[tri]
        poly = Poly3DCollection([verts], alpha=1.0, facecolor=facecolor,
                                edgecolor=edgecolor, linewidth=0.5)
        ax.add_collection3d(poly)

    pad = 0.1
    for dim in range(3):
        lo, hi = proj[:, dim].min(), proj[:, dim].max()
        margin = (hi - lo) * pad
        if dim == 0: ax.set_xlim(lo - margin, hi + margin)
        elif dim == 1: ax.set_ylim(lo - margin, hi + margin)
        else: ax.set_zlim(lo - margin, hi + margin)

    ax.set_title(title, fontsize=11)
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')


def main():
    # ── 10→3 ──
    gens10, ico, faces = build_10_3_generators()

    print("=" * 70)
    print("10→3 GENERATORS: Icosahedral three-fold axes")
    print("  Rhombic enneacontahedron, I_h symmetry (order 120)")
    print("  f-vector: (92, 180, 90)")
    print("  Free parameters: 0")
    print("=" * 70)
    print()
    print("Icosahedron vertices: cyclic permutations of (0, ±1, ±φ)")
    print(f"where φ = (1+√5)/2 = {PHI:.10f}")
    print()
    print("Generators (unit vectors along face centers, one per antipodal pair):")
    for i, g in enumerate(gens10):
        print(f"  g{i:2d} = ({g[0]:+.8f}, {g[1]:+.8f}, {g[2]:+.8f})")

    # Verify tight frame
    P10 = gens10.T
    G10 = P10 @ P10.T
    print(f"\nPP^T = {G10[0,0]:.6f} · I  (should be 10/3 = {10/3:.6f})")

    vol10 = sum(abs(det(gens10[list(idx)])) for idx in itertools.combinations(range(10), 3))
    print(f"Raw volume = {vol10:.10f}")
    print(f"Normalized = {vol10 / (10/3)**1.5:.10f}")

    # ── 8→3 ──
    gens8, a, b = build_8_3_generators()

    print()
    print("=" * 70)
    print("8→3 GENERATORS: 7 vectors on cone + 1 axial")
    print("  D_7 × Z_2 symmetry (order 28)")
    print("  f-vector: (58, 112, 56)")
    print("  Free parameters: 1 (θ ≈ 1.06578)")
    print("=" * 70)
    print()
    print("Symbolic parameterization:")
    print("  g_k = a·(sinθ·cos(2πk/7), sinθ·sin(2πk/7), cosθ)  for k=0..6")
    print("  g_7 = (0, 0, b)")
    print()
    print("  a² = 16 / (21·sin²θ)")
    print("  b² = 8·(1 - 3·cos²θ) / (3·sin²θ)")
    print()
    print(f"  θ  = {1.065776577865147:.15f}")
    print(f"  a  = {a:.15f}")
    print(f"  b  = {b:.15f}")
    print()
    print("Generators:")
    for i, g in enumerate(gens8):
        label = f"k={i}" if i < 7 else "axial"
        print(f"  g{i} = ({g[0]:+.10f}, {g[1]:+.10f}, {g[2]:+.10f})  [{label}]")

    vol8 = sum(abs(det(gens8[list(idx)])) for idx in itertools.combinations(range(8), 3))
    print(f"\nRaw volume = {vol8:.10f}")
    print(f"Normalized = {vol8 / (8/3)**1.5:.10f}")

    # ── Plots ──
    proj10 = build_zonotope(gens10, 10)
    proj8 = build_zonotope(gens8, 8)

    # 10→3: Rhombic enneacontahedron
    fig1 = plt.figure(figsize=(8, 7))
    ax1 = fig1.add_subplot(111, projection='3d')
    plot_zonotope(ax1, proj10,
                  "10-cube → 3D: Rhombic enneacontahedron\n"
                  "I_h symmetry (order 120), 90 rhombic faces, f = (92, 180, 90)",
                  'steelblue', 'navy')
    ax1.view_init(elev=20, azim=35)
    plt.tight_layout()
    plt.savefig('zonotope_10_3.png', dpi=150, bbox_inches='tight')
    print("\nSaved: zonotope_10_3.png")
    plt.close()

    # 10→3: top-down view along z-axis (5-fold symmetry visible)
    fig1b = plt.figure(figsize=(8, 7))
    ax1b = fig1b.add_subplot(111, projection='3d')
    plot_zonotope(ax1b, proj10,
                  "10-cube → 3D: Rhombic enneacontahedron (top view)\n"
                  "I_h symmetry (order 120)",
                  'steelblue', 'navy')
    ax1b.view_init(elev=90, azim=0)
    plt.tight_layout()
    plt.savefig('zonotope_10_3_top.png', dpi=150, bbox_inches='tight')
    print("Saved: zonotope_10_3_top.png")
    plt.close()

    # 8→3: D_7 x Z_2 zonotope
    fig2 = plt.figure(figsize=(8, 7))
    ax2 = fig2.add_subplot(111, projection='3d')
    plot_zonotope(ax2, proj8,
                  "8-cube → 3D: D₇×Z₂ zonotope\n"
                  "Order 28, 56 parallelogram faces, f = (58, 112, 56)",
                  'coral', 'darkred')
    ax2.view_init(elev=20, azim=35)
    plt.tight_layout()
    plt.savefig('zonotope_8_3.png', dpi=150, bbox_inches='tight')
    print("Saved: zonotope_8_3.png")
    plt.close()

    # 8→3: top-down view along 7-fold axis
    fig3 = plt.figure(figsize=(8, 7))
    ax3 = fig3.add_subplot(111, projection='3d')
    plot_zonotope(ax3, proj8,
                  "8-cube → 3D: view along 7-fold axis (z)\n"
                  "D₇×Z₂ symmetry (order 28)",
                  'coral', 'darkred')
    ax3.view_init(elev=90, azim=0)
    plt.tight_layout()
    plt.savefig('zonotope_8_3_top.png', dpi=150, bbox_inches='tight')
    print("Saved: zonotope_8_3_top.png")


if __name__ == '__main__':
    main()
