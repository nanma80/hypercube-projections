"""Draw the max-area 2D projection of the 4_21 polytope.

Shows the regular hexagon hull (6 vertices at radius sqrt(2)) plus all 240
projected E8 roots, color-coded by multiplicity. Saves to PNG.
"""
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon, Circle
from matplotlib.collections import PatchCollection
from pathlib import Path
from math import sqrt

# ---- Build E8 roots ----
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si in (-1, 1):
        for sj in (-1, 1):
            v = [0.0]*8; v[i] = si; v[j] = sj; roots.append(v)
for s in itertools.product((-1, 1), repeat=8):
    if int(np.prod(s)) == 1:
        roots.append([0.5 * x for x in s])
ROOTS = np.array(roots)

# ---- A2 subsystem plane: r1 = e1+e2, r2 = e1-e3 ----
r1 = np.array([1, 1, 0, 0, 0, 0, 0, 0], dtype=float)
r2 = np.array([1, 0, -1, 0, 0, 0, 0, 0], dtype=float)
# Orthonormal basis
u1 = r1 / np.linalg.norm(r1)
u2 = r2 - (r2 @ u1) * u1
u2 = u2 / np.linalg.norm(u2)
Q = np.column_stack([u1, u2])  # 8x2

pts = ROOTS @ Q  # 240 x 2

# ---- Group by distinct 2D points and count multiplicities ----
pts_round = np.round(pts, 6)
point_map = {}
for p in pts_round:
    key = tuple(p)
    point_map[key] = point_map.get(key, 0) + 1

distinct_pts = np.array(list(point_map.keys()))
mults = np.array(list(point_map.values()))

print(f"240 roots -> {len(distinct_pts)} distinct 2D points")
print(f"Multiplicity distribution:")
for m in sorted(set(mults)):
    n = sum(mults == m)
    print(f"  mult {m}: {n} points")

# ---- Identify hull vertices ----
hull_mask = np.array([abs(np.linalg.norm(p) - sqrt(2)) < 1e-4 for p in distinct_pts])

# ---- Plot ----
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.set_aspect('equal')

# Draw hexagon outline
hex_verts = []
for k in range(6):
    angle = k * np.pi / 3 + np.pi / 6  # start from 30° to align with our vertices
    hex_verts.append([sqrt(2) * np.cos(angle), sqrt(2) * np.sin(angle)])
hex_verts.append(hex_verts[0])  # close
hx, hy = zip(*hex_verts)

# Find actual hull vertex angles to draw the hexagon correctly
hull_pts = distinct_pts[hull_mask]
hull_angles = sorted(np.arctan2(p[1], p[0]) for p in hull_pts)
hex_outline = []
for a in hull_angles:
    hex_outline.append([sqrt(2) * np.cos(a), sqrt(2) * np.sin(a)])
hex_outline.append(hex_outline[0])
hox, hoy = zip(*hex_outline)
ax.plot(hox, hoy, 'k-', linewidth=1.5, zorder=2)

# Fill hexagon with very light color
from matplotlib.patches import Polygon as MplPolygon
hex_poly = MplPolygon(hex_outline[:-1], closed=True, facecolor='#f0f4ff',
                       edgecolor='none', zorder=1)
ax.add_patch(hex_poly)

# ---- Draw projected E8 edges ----
# E8 edges: pairs of roots with inner product 1 (nearest neighbors)
# Project all edges to 2D and draw, grouping by multiplicity
from collections import Counter
edge_2d_count = Counter()
for i in range(len(ROOTS)):
    for j in range(i + 1, len(ROOTS)):
        if abs(ROOTS[i] @ ROOTS[j] - 1.0) < 1e-8:
            p_i = tuple(np.round(pts[i], 6))
            p_j = tuple(np.round(pts[j], 6))
            if p_i != p_j:  # skip zero-length edges
                key = tuple(sorted([p_i, p_j]))
                edge_2d_count[key] += 1

print(f"Distinct projected edges: {len(edge_2d_count)}")
mult_vals_e = sorted(set(edge_2d_count.values()))
print(f"Edge multiplicity values: {mult_vals_e}")
for m in mult_vals_e:
    n = sum(1 for v in edge_2d_count.values() if v == m)
    print(f"  mult {m}: {n} edges")

# Draw edges — scale linewidth by log of multiplicity for visibility
max_emult = max(edge_2d_count.values())
for (p1, p2), emult in sorted(edge_2d_count.items(), key=lambda x: x[1]):
    t = np.log1p(emult) / np.log1p(max_emult)
    lw = 1.0 + 2.5 * t
    alpha = 0.4 + 0.5 * t
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
            color='#4477aa', linewidth=lw, alpha=alpha, zorder=1.5)

# Color map by multiplicity
max_mult = max(mults)
cmap = plt.cm.YlOrRd

# Sort by multiplicity (draw low mult first, high on top)
order = np.argsort(mults)

# Draw interior points first, hull points last
for idx in order:
    p = distinct_pts[idx]
    m = mults[idx]
    r = np.linalg.norm(p)

    if r < 1e-4:
        # Origin: special marker
        size = 120
        color = '#333333'
        marker = 'D'
        zorder = 5
    elif abs(r - sqrt(2)) < 1e-4:
        # Hull vertex
        size = 100
        color = '#d62728'
        marker = 'o'
        zorder = 6
    else:
        # Interior point
        t = m / max_mult
        size = 30 + 50 * t
        color = cmap(0.2 + 0.6 * t)
        marker = 'o'
        zorder = 3

    ax.scatter(p[0], p[1], s=size, c=[color], marker=marker,
               edgecolors='black', linewidths=0.5, zorder=zorder)

# Label hull vertices with their E8 root
for p in hull_pts:
    angle = np.arctan2(p[1], p[0])
    offset = 0.18
    lx = p[0] + offset * np.cos(angle)
    ly = p[1] + offset * np.sin(angle)

    # Find the E8 root
    for ir in range(len(ROOTS)):
        proj = ROOTS[ir] @ Q
        if np.allclose(proj, p, atol=1e-4):
            r8 = ROOTS[ir]
            # Format nicely
            nonzero = [(k, r8[k]) for k in range(8) if abs(r8[k]) > 0.01]
            if all(abs(v) == 1.0 for _, v in nonzero):
                label = "".join(f"{'+'if v>0 else '−'}e{k+1}" for k, v in nonzero)
            else:
                label = ""
            if label:
                ax.annotate(label, (p[0], p[1]), (lx, ly),
                           fontsize=7, ha='center', va='center',
                           color='#666666')
            break

# Annotations
ax.annotate(f"origin\n(72 roots)", (0, 0), (0.35, -0.35),
            fontsize=8, ha='left', color='#333333',
            arrowprops=dict(arrowstyle='->', color='#999999', lw=0.8))

# Add multiplicity legend with matching colors/markers
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728',
           markeredgecolor='black', markersize=9, label='hull vertices (mult 1, 6 pts)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(0.2 + 0.6 * 27/max_mult),
           markeredgecolor='black', markersize=7, label='interior (mult 27, 6 pts)'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#333333',
           markeredgecolor='black', markersize=8, label='origin (mult 72, 1 pt)'),
]
ax.legend(handles=legend_handles, loc='upper left', fontsize=8, framealpha=0.9)

# Title and labels
ax.set_title("Max-area 2D projection of 4₂₁ polytope (E₈ roots)\n"
             f"Regular hexagon on A₂ plane, area = 3√3 ≈ {3*sqrt(3):.4f}",
             fontsize=11, pad=12)

ax.set_xlim(-1.9, 1.9)
ax.set_ylim(-1.9, 1.9)
ax.set_axis_off()

plt.tight_layout()
out = Path(__file__).parent / "max_area_421_to_2d.png"
plt.savefig(out, dpi=200, bbox_inches='tight')
print(f"Saved {out}")
plt.close()
