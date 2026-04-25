"""Generate STL of all distinct 4_21 struts projected to 3D.

Builds: 49 vertex balls + 432 edge cylinders.
Ball radius takes one of 3 discrete sizes (1x, 2x, 3x base) by ranking
the distinct multiplicity values: smallest mult -> 1x, middle -> 2x,
largest -> 3x.
"""
import numpy as np
import itertools
import trimesh
from scipy.spatial import ConvexHull
from collections import Counter
from pathlib import Path

s2 = np.sqrt(2.0)
s23 = np.sqrt(2.0/3.0)
s21 = np.sqrt(21.0)
s6 = np.sqrt(6.0)

# Closed-form 8 generators
gens = np.array([
    [+1/s2, 0,      3/(4*s21)],
    [-1/s2, 0,      3/(4*s21)],
    [0,    -1/s6,   7/(4*s21)],
    [0,    +1/s6, -13/(4*s21)],
    [0,    +1/s6,   5/(4*s21)],
    [0,    +1/s6,   5/(4*s21)],
    [0,    +1/s6,   5/(4*s21)],
    [0,    +1/s6,   5/(4*s21)],
])

# 240 E_8 roots
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product([-1, 1], repeat=2):
        r = [0.0]*8; r[i] = si; r[j] = sj; roots.append(r)
for s in itertools.product([-1, 1], repeat=8):
    if int(np.prod(s)) == 1:
        roots.append([0.5*x for x in s])
roots = np.array(roots)

# 4_21 edges in 8D
G = roots @ roots.T
edges_8d = [(i, j) for i in range(240) for j in range(i+1, 240) if abs(G[i, j] - 1) < 1e-9]
print(f"E_8 roots: {len(roots)},  4_21 edges in 8D: {len(edges_8d)}")

# Project and dedupe
P = roots @ gens
unique_pts = []
pt_idx = []
for p in P:
    found = -1
    for k, u in enumerate(unique_pts):
        if np.linalg.norm(p - u) < 1e-8:
            found = k; break
    if found < 0:
        unique_pts.append(p); found = len(unique_pts) - 1
    pt_idx.append(found)
unique_pts = np.array(unique_pts)
mult = Counter(pt_idx)
print(f"Distinct projected points: {len(unique_pts)}")
print(f"Multiplicity distribution: {dict(Counter(mult.values()))}")

# Distinct projected nonzero segments + multiplicity
seg_mult = Counter()
for i, j in edges_8d:
    a, b = pt_idx[i], pt_idx[j]
    if a == b: continue
    seg_mult[tuple(sorted([a, b]))] += 1
print(f"Distinct projected struts: {len(seg_mult)}")
print(f"Strut multiplicity distribution: {dict(Counter(seg_mult.values()))}")

hull = ConvexHull(unique_pts)
print(f"Hull vertices: {len(hull.vertices)}, hull volume: {hull.volume:.6f}")

# Build mesh; ball radius = base * rank-of-multiplicity (1x, 2x, 3x, ...)
# Strut radius = base * (1.0 + 0.5 * (rank - 1)) = base * 0.5 * (rank + 1).
base_ball_r  = 0.07
base_stick_r = 0.012

sorted_pt_mults = sorted(set(mult.values()))
ptmult_to_size  = {m: 1.0 + 0.5 * i for i, m in enumerate(sorted_pt_mults)}
sorted_seg_mults = sorted(set(seg_mult.values()))
segmult_to_size  = {m: 1.0 + 0.5 * i for i, m in enumerate(sorted_seg_mults)}

meshes = []
for k, v in enumerate(unique_pts):
    r = base_ball_r * ptmult_to_size[mult[k]]
    s_ = trimesh.creation.icosphere(subdivisions=2, radius=r)
    s_.apply_translation(v); meshes.append(s_)

z_axis = np.array([0., 0., 1.])
for (a, b), m in seg_mult.items():
    p, q = unique_pts[a], unique_pts[b]
    seg = q - p
    L = np.linalg.norm(seg)
    r_stick = base_stick_r * segmult_to_size[m]
    cyl = trimesh.creation.cylinder(radius=r_stick, height=L, sections=12)
    direction = seg / L
    if np.allclose(direction, z_axis):
        R = np.eye(4)
    elif np.allclose(direction, -z_axis):
        R = trimesh.transformations.rotation_matrix(np.pi, [1, 0, 0])
    else:
        axis = np.cross(z_axis, direction); axis /= np.linalg.norm(axis)
        angle = np.arccos(np.clip(np.dot(z_axis, direction), -1, 1))
        R = trimesh.transformations.rotation_matrix(angle, axis)
    cyl.apply_transform(R)
    cyl.apply_translation((p + q) / 2.0)
    meshes.append(cyl)

mesh = trimesh.util.concatenate(meshes)
out = Path(__file__).parent / "max_volume_421_to_3d_all_struts.stl"
mesh.export(out)
print(f"\nWrote {out}")
print(f"  triangles: {len(mesh.faces)}")
print(f"  balls: {len(unique_pts)} with radii base * (1.0, 1.5, ...), base = {base_ball_r}")
for m in sorted_pt_mults:
    sz = ptmult_to_size[m]
    n = sum(1 for k in range(len(unique_pts)) if mult[k] == m)
    print(f"    pt mult={m:2d} -> {sz}x: {n:3d} balls, r = {base_ball_r * sz:.4f}")
print(f"  cylinders: {len(seg_mult)} with radii base * (1.0, 1.5, ...), base = {base_stick_r}")
for m in sorted_seg_mults:
    sz = segmult_to_size[m]
    n = sum(1 for v in seg_mult.values() if v == m)
    print(f"    seg mult={m:3d} -> {sz}x: {n:3d} cylinders, r = {base_stick_r * sz:.4f}")

