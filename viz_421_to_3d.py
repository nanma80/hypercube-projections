#!/usr/bin/env python3
"""
Visualize the max-volume 3D shadow of the 4_21 polytope and dump the
8 generator vectors that define it.

The 240 4_21 vertices in 3D are linear combinations of 8 generator
vectors g_0..g_7 (the projections of the 8 standard basis vectors of R^8):
  - integer roots:    +-g_i +- g_j        (i!=j)         -> 112 vertices
  - half-integer:     (1/2) sum_k s_k g_k (even-parity s) -> 128 vertices

After basin-hopping, we canonicalize the orientation so:
  - The principal C_3 axis is aligned with the +z axis
  - One C_2 axis is aligned with the +x axis (so a sigma_d mirror is the xz plane)

Outputs:
  - max_volume_421_to_3d_generators.txt : the 8 generators
  - displays a 3D matplotlib plot
"""
import argparse
import itertools
import numpy as np
from scipy.linalg import orth
from scipy.spatial import ConvexHull
from scipy.optimize import basinhopping


def get_4_21_roots():
    roots = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (-1, 1):
            for sj in (-1, 1):
                v = [0.0] * 8
                v[i] = si
                v[j] = sj
                roots.append(v)
    for s in itertools.product((-1, 1), repeat=8):
        if int(np.prod(s)) == 1:
            roots.append([0.5 * x for x in s])
    return np.array(roots)


ROOTS = get_4_21_roots()


def neg_vol(p):
    Q = orth(p.reshape(3, 8).T)
    if Q.shape[1] < 3:
        return 0.0
    try:
        return -ConvexHull(ROOTS @ Q).volume
    except Exception:
        return 0.0


def find_optimum(n_trials=20, niter=80, seed=42):
    rng = np.random.default_rng(seed)
    best, bestM = 0.0, None
    for t in range(n_trials):
        x0 = rng.standard_normal(24)
        res = basinhopping(
            neg_vol, x0, niter=niter,
            minimizer_kwargs={'method': 'L-BFGS-B',
                              'options': {'eps': 1e-4}},
            seed=int(rng.integers(1 << 30)), disp=False,
        )
        v = -res.fun
        if v > best:
            best, bestM = v, res.x.reshape(3, 8)
            print(f"  trial {t+1:3d}: new best = {v:.10f}")
    return best, bestM


def find_symmetry_group(verts, atol=1e-4):
    """Return list of 3x3 orthogonal matrices preserving verts (set)."""
    n = len(verts)
    # Find a non-degenerate base triple
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                M = verts[[i, j, k]]
                if abs(np.linalg.det(M)) > 0.1:
                    base_idx, base_M = (i, j, k), M
                    break
            else:
                continue
            break
        else:
            continue
        break
    G_base = base_M @ base_M.T
    inv_base = np.linalg.inv(base_M)

    sym = []
    for a in range(n):
        for b in range(n):
            if b == a:
                continue
            for c in range(n):
                if c in (a, b):
                    continue
                M2 = verts[[a, b, c]]
                if not np.allclose(M2 @ M2.T, G_base, atol=atol):
                    continue
                R = inv_base @ M2
                if not np.allclose(R.T @ R, np.eye(3), atol=atol):
                    continue
                ok = True
                Vt = verts @ R
                for p in Vt:
                    if not any(np.allclose(p, q, atol=atol) for q in verts):
                        ok = False
                        break
                if ok:
                    sym.append(R)
    return sym


def canonicalize(Q):
    """Rotate Q (8x3) so the C_3 axis is +z and a C_2 axis is +x.

    Returns Q_canon (8x3).
    """
    pts = ROOTS @ Q
    hull = ConvexHull(pts)
    HV = pts[hull.vertices]
    sym = find_symmetry_group(HV, atol=1e-3)

    # Find C_3 rotation axis
    c3_axis = None
    for R in sym:
        if abs(np.linalg.det(R) - 1) < 1e-2:
            tr = np.trace(R)
            ang = np.arccos(np.clip((tr - 1) / 2, -1, 1))
            if abs(ang - 2 * np.pi / 3) < 0.05:
                # axis = eigenvector with eigenvalue 1
                w, V = np.linalg.eig(R)
                for i, e in enumerate(w):
                    if abs(e - 1) < 0.05:
                        ax = np.real(V[:, i])
                        ax = ax / np.linalg.norm(ax)
                        c3_axis = ax
                        break
                if c3_axis is not None:
                    break

    if c3_axis is None:
        print("  WARNING: no C_3 axis found; skipping canonicalization")
        return Q

    # Rotate so c3_axis -> +z
    z = np.array([0., 0., 1.])
    v = np.cross(c3_axis, z)
    s = np.linalg.norm(v)
    c = c3_axis @ z
    if s < 1e-9:
        R1 = np.eye(3) if c > 0 else np.diag([1., -1., -1.])
    else:
        K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R1 = np.eye(3) + K + K @ K * ((1 - c) / (s * s))

    Q1 = Q @ R1.T

    # Recompute symmetry group of rotated Q1's vertex set
    pts1 = ROOTS @ Q1
    hull1 = ConvexHull(pts1)
    HV1 = pts1[hull1.vertices]
    sym1 = find_symmetry_group(HV1, atol=1e-3)

    # Find a C_2 axis perpendicular to z
    c2_axis = None
    for R in sym1:
        if abs(np.linalg.det(R) - 1) < 1e-2:
            tr = np.trace(R)
            ang = np.arccos(np.clip((tr - 1) / 2, -1, 1))
            if abs(ang - np.pi) < 0.05:
                w, V = np.linalg.eig(R)
                for i, e in enumerate(w):
                    if abs(e - 1) < 0.05:
                        ax = np.real(V[:, i])
                        ax = ax / np.linalg.norm(ax)
                        if abs(ax[2]) < 0.05:  # perpendicular to z
                            c2_axis = ax
                            break
                if c2_axis is not None:
                    break

    if c2_axis is None:
        return Q1

    # Rotate around z so c2_axis -> +x
    theta = -np.arctan2(c2_axis[1], c2_axis[0])
    cz, sz = np.cos(theta), np.sin(theta)
    R2 = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Q1 @ R2.T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=15)
    ap.add_argument('--niter', type=int, default=80)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--save', default='AI_assisted_analysis/4_21_to_3d/'
                                       'max_volume_421_to_3d_generators.txt')
    ap.add_argument('--no-plot', action='store_true')
    args = ap.parse_args()

    print(f"Optimizing {args.trials} trials...")
    vol, M = find_optimum(args.trials, args.niter, args.seed)
    Q = orth(M.T)  # 8x3 orthonormal cols
    print(f"\nMax volume = {vol:.10f}")

    print("\nCanonicalizing orientation (C_3 axis -> +z, C_2 -> +x)...")
    Q = canonicalize(Q)
    pts = ROOTS @ Q
    hull = ConvexHull(pts)
    print(f"After canonicalization: vol = {hull.volume:.10f}, |V|={len(hull.vertices)}")

    # The 8 generators are the rows of Q
    generators = Q
    print("\n8 generator vectors (rows of canonicalized Q):")
    for i, g in enumerate(generators):
        print(f"  g_{i} = ({g[0]:+.10f}, {g[1]:+.10f}, {g[2]:+.10f})")

    # Look at structure: which generators are equal?
    print("\nGenerator equality classes (4 axes coalesce):")
    classes = []
    used = [False] * 8
    for i in range(8):
        if used[i]:
            continue
        cls = [i]
        used[i] = True
        for j in range(i+1, 8):
            if not used[j] and np.allclose(generators[i], generators[j], atol=1e-4):
                cls.append(j)
                used[j] = True
        classes.append(cls)
    for c in classes:
        gi = generators[c[0]]
        print(f"  axes {c}: g = ({gi[0]:+.6f}, {gi[1]:+.6f}, {gi[2]:+.6f}), "
              f"|g| = {np.linalg.norm(gi):.6f}")

    # Save generators
    import os
    os.makedirs(os.path.dirname(args.save), exist_ok=True)
    with open(args.save, 'w') as f:
        f.write("# Max-volume 3D projection of the 4_21 polytope\n")
        f.write(f"# Volume: {vol:.10f}\n")
        f.write("# Orthonormal projection basis Q (8x3): rows are the 8 generators g_0..g_7\n")
        f.write("# The 240 4_21 vertices in 3D are:\n")
        f.write("#   integer roots:    +-g_i +- g_j  (i != j)              -> 112\n")
        f.write("#   half-integer:     (1/2) sum_k s_k g_k (even parity)   -> 128\n")
        f.write("# Canonicalized: C_3 axis = +z, one C_2 axis = +x\n")
        f.write("# Symmetry: D_3d (order 12)\n")
        f.write(f"# Hull: {len(hull.vertices)} vertices, {len(hull.simplices)} triangulated facets\n")
        f.write("\n")
        for i, g in enumerate(generators):
            f.write(f"g_{i}: {g[0]:+.15f}, {g[1]:+.15f}, {g[2]:+.15f}\n")
        f.write("\n# Hull vertex coordinates (after canonicalization):\n")
        HV = pts[hull.vertices]
        # sort by z then by angle in xy
        order = sorted(range(len(HV)),
                       key=lambda k: (round(HV[k][2], 4),
                                      np.arctan2(HV[k][1], HV[k][0])))
        for idx, k in enumerate(order):
            v = HV[k]
            f.write(f"v_{idx:02d}: {v[0]:+.10f}, {v[1]:+.10f}, {v[2]:+.10f}  "
                    f"|v| = {np.linalg.norm(v):.6f}\n")
    print(f"\nSaved generators + hull vertices to: {args.save}")

    if args.no_plot:
        return

    # Visualization
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    from collections import Counter

    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection='3d')

    # Group hull vertices by norm to color-code orbits
    HV = pts[hull.vertices]
    norms = np.linalg.norm(HV, axis=1)
    unique_norms = sorted({float(round(n, 4)) for n in norms})
    colors = {n: c for n, c in zip(unique_norms, ['red', 'orange', 'royalblue'])}
    for v, n in zip(HV, norms):
        ax.scatter(*v, color=colors[float(round(n, 4))], s=70,
                   edgecolor='black', linewidth=0.5, zorder=5)

    # Draw faces (triangulated by ConvexHull) -- merge coplanar triangles for cleaner look
    # Collect triangles with face color by their normal direction class
    tris = []
    for s in hull.simplices:
        tri = [pts[s[0]], pts[s[1]], pts[s[2]]]
        tris.append(tri)
    poly = Poly3DCollection(tris, alpha=0.18, facecolor='lightsteelblue',
                            edgecolor='steelblue', linewidth=0.4)
    ax.add_collection3d(poly)

    # Draw real edges (those bordering >= 2 distinct face planes)
    eqs = hull.equations
    plane_keys = []
    for n_, d_ in zip(eqs[:, :3], eqs[:, 3]):
        if (n_[0] < -1e-9
            or (abs(n_[0]) < 1e-9 and n_[1] < -1e-9)
            or (abs(n_[0]) < 1e-9 and abs(n_[1]) < 1e-9 and n_[2] < 0)):
            n_, d_ = -n_, -d_
        plane_keys.append((tuple(np.round(n_, 5)), round(d_, 5)))
    plane_to_facets = {}
    for i, k in enumerate(plane_keys):
        plane_to_facets.setdefault(k, []).append(i)

    real_edges = set()
    face_vertex_sets = []
    for k, idxs in plane_to_facets.items():
        vs = set()
        for i in idxs:
            vs.update(int(x) for x in hull.simplices[i])
        face_vertex_sets.append(vs)
    for fvs in face_vertex_sets:
        for a, b in itertools.combinations(fvs, 2):
            cnt = sum(1 for fs in face_vertex_sets if a in fs and b in fs)
            if cnt >= 2:
                real_edges.add((min(a, b), max(a, b)))

    elens_raw = [np.linalg.norm(pts[a] - pts[b]) for a, b in real_edges]
    elens = [float(round(l, 3)) for l in elens_raw]
    elen_classes = sorted(set(elens))
    palette = ['black', 'green', 'purple', 'crimson', 'gold', 'teal']
    edge_colors = {l: palette[i % len(palette)] for i, l in enumerate(elen_classes)}
    segs = [[pts[a], pts[b]] for a, b in real_edges]
    seg_colors = [edge_colors[el] for el in elens]
    lc = Line3DCollection(segs, colors=seg_colors, linewidth=1.6)
    ax.add_collection3d(lc)

    # Legend
    from matplotlib.lines import Line2D
    legend_elems = []
    norm_counts = Counter(float(round(n, 4)) for n in norms)
    for n, c in colors.items():
        legend_elems.append(Line2D([0], [0], marker='o', color='w',
                                   markerfacecolor=c, markeredgecolor='k',
                                   markersize=10,
                                   label=f"|v|={n:.4f} ({norm_counts[n]} verts)"))
    elen_counts = Counter(elens)
    for l, c in edge_colors.items():
        legend_elems.append(Line2D([0], [0], color=c, linewidth=2.5,
                                   label=f"edge len={l:.4f} ({elen_counts[l]} edges)"))
    ax.legend(handles=legend_elems, loc='upper left', fontsize=8,
              framealpha=0.9)

    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_title(f'Max-volume 4_21 -> 3D shadow\n'
                 f'18 vertices, 48 edges, 32 triangular faces, '
                 f'D_3d symmetry (vol={vol:.6f})')

    # Equal aspect
    extent = max(np.ptp(pts[:, k]) for k in range(3)) / 2
    ctr = pts.mean(axis=0)
    for f, c in zip([ax.set_xlim, ax.set_ylim, ax.set_zlim], ctr):
        f(c - extent, c + extent)
    try:
        ax.set_box_aspect([1, 1, 1])
    except Exception:
        pass

    plt.tight_layout()
    out_png = args.save.replace('.txt', '.png')
    plt.savefig(out_png, dpi=140)
    print(f"Saved plot to: {out_png}")
    plt.show()


if __name__ == '__main__':
    main()
