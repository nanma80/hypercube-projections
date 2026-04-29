"""Generate STL wireframes of the max-volume 4_21->4D hull (F4 root system).

The 48 hull vertices are the F4 root system with closed-form coordinates:

  Orbit A — 24 "long" roots, |v|^2 = (5+sqrt(7))/4:
    All permutations of (+-a, +-a, 0, 0) where a = sqrt((5+sqrt(7))/8)

  Orbit B — 24 "short" roots, |v|^2 = (3+sqrt(7))/4:
    All permutations of (+-b, 0, 0, 0)           (8 vertices)
    All sign choices (+-b/2, +-b/2, +-b/2, +-b/2) (16 vertices)
    where b = sqrt((3+sqrt(7))/4)

Two orthographic 4D->3D projections are generated, one for each type of
root at the "pole" (i.e. along the projection direction):

  Short-root view: project along (1,1,1,1)/2 (an orbit B direction).
    - An orbit B vertex sits at the pole/antipole.
    - Preserved symmetry: S4 (all permutations of coords), order 24.
    - 3D complement basis: (1,-1,0,0)/sqrt(2), (1,1,-2,0)/sqrt(6),
      (1,1,1,-3)/sqrt(12).

  Long-root view: project along (1,1,0,0)/sqrt(2) (an orbit A direction).
    - An orbit A vertex sits at the pole/antipole.
    - Preserved symmetry: D4 x Z2 (reflections in all 3 axes + yz swap
      + C4 in yz plane + inversion), order 16.
    - 3D complement basis: (1,-1,0,0)/sqrt(2), (0,0,1,0), (0,0,0,1).

Both views have 7 depth layers (symmetric about 0) and produce 33 distinct
3D points from the 48 vertices. Orbit A balls are larger (0.06) and orbit B
balls smaller (0.045).

Outputs:
  max_volume_421_to_4d_short_root_view.stl
  max_volume_421_to_4d_long_root_view.stl
"""
import itertools
import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from pathlib import Path

# ---- Build the 48 F4 vertices in closed form ----
s7 = np.sqrt(7.0)
a = np.sqrt((5 + s7) / 8.0)   # for orbit A: perms of (+-a, +-a, 0, 0)
b = np.sqrt((3 + s7) / 4.0)   # for orbit B: perms of (+-b, 0, 0, 0) + (+-b/2)^4

orbit_A = []
for i, j in itertools.combinations(range(4), 2):
    for si in (-1, 1):
        for sj in (-1, 1):
            v = [0.0]*4
            v[i] = si * a
            v[j] = sj * a
            orbit_A.append(v)

orbit_B = []
# axis vertices
for i in range(4):
    for s in (-1, 1):
        v = [0.0]*4
        v[i] = s * b
        orbit_B.append(v)
# all-half vertices
for signs in itertools.product((-1, 1), repeat=4):
    orbit_B.append([s * b / 2 for s in signs])

orbit_A = np.array(orbit_A)
orbit_B = np.array(orbit_B)
V4 = np.vstack([orbit_A, orbit_B])
print(f"4D vertices: {len(V4)} (24 + 24)")

# Verify hull
hull4 = ConvexHull(V4)
V_target = (17 + 7*s7) / 4
print(f"4D hull volume: {hull4.volume:.10f}")
print(f"Target (17+7*sqrt(7))/4: {V_target:.10f}")
print(f"Hull has {len(hull4.vertices)} vertices, {len(hull4.simplices)} simplices")

# ---- Get the 4D edges (ridges shared by exactly 2 facets) ----
# In 4D, edges are 1-faces. ConvexHull.simplices are tetrahedra (3-simplices).
# An edge is a pair of vertices appearing in at least 2 simplices.
# More precisely: edges of the polytope are pairs of vertices that are
# connected by a 1-face of the convex hull.
edge_count = {}
for simplex in hull4.simplices:
    for i, j in itertools.combinations(simplex, 2):
        key = (min(i, j), max(i, j))
        edge_count[key] = edge_count.get(key, 0) + 1

# An edge of a 4-polytope appears in at least 2 tetrahedral facets
edges_4d = [e for e, c in edge_count.items() if c >= 2]
print(f"4D polytope edges: {len(edges_4d)}")

# Classify edges
n_AA = sum(1 for i, j in edges_4d if i < 24 and j < 24)
n_BB = sum(1 for i, j in edges_4d if i >= 24 and j >= 24)
n_AB = sum(1 for i, j in edges_4d if (i < 24) != (j < 24))
print(f"  A-A edges: {n_AA}, B-B edges: {n_BB}, A-B edges: {n_AB}")

# Edge lengths
elens = {}
for i, j in edges_4d:
    d = round(np.linalg.norm(V4[i] - V4[j]), 4)
    elens[d] = elens.get(d, 0) + 1
print(f"  Edge length distribution: {dict(sorted(elens.items()))}")

# ---- Define both projection views ----
views = [
    {
        "name": "short_root_view",
        "desc": "Short-root pole: (1,1,1,1)/2 - S4 symmetry, orbit B at pole",
        "basis": np.column_stack([
            np.array([1, -1, 0, 0]) / np.sqrt(2),
            np.array([1, 1, -2, 0]) / np.sqrt(6),
            np.array([1, 1, 1, -3]) / np.sqrt(12),
        ]),
    },
    {
        "name": "long_root_view",
        "desc": "Long-root pole: (1,1,0,0)/sqrt(2) - D4xZ2 symmetry, orbit A at pole",
        "basis": np.column_stack([
            np.array([1, -1, 0, 0]) / np.sqrt(2),
            np.array([0, 0, 1, 0]),
            np.array([0, 0, 0, 1]),
        ]),
    },
]

# ---- Build STL for each view ----
ball_radius_A = 0.06   # larger for orbit A
ball_radius_B = 0.045  # smaller for orbit B
stick_radius = 0.015
z_axis = np.array([0., 0., 1.])
out_dir = Path(__file__).parent

for view in views:
    print(f"\n--- {view['desc']} ---")
    V3 = V4 @ view["basis"]
    print(f"3D projected vertices: {len(V3)}")

    meshes = []
    for idx in range(len(V3)):
        r = ball_radius_A if idx < 24 else ball_radius_B
        s = trimesh.creation.icosphere(subdivisions=2, radius=r)
        s.apply_translation(V3[idx])
        meshes.append(s)

    for i, j in edges_4d:
        p, q = V3[i], V3[j]
        seg = q - p
        L = np.linalg.norm(seg)
        if L < 1e-6:
            continue
        cyl = trimesh.creation.cylinder(radius=stick_radius, height=L, sections=12)
        direction = seg / L
        if np.allclose(direction, z_axis):
            R = np.eye(4)
        elif np.allclose(direction, -z_axis):
            R = trimesh.transformations.rotation_matrix(np.pi, [1, 0, 0])
        else:
            axis = np.cross(z_axis, direction)
            axis /= np.linalg.norm(axis)
            angle = np.arccos(np.clip(np.dot(z_axis, direction), -1, 1))
            R = trimesh.transformations.rotation_matrix(angle, axis)
        cyl.apply_transform(R)
        cyl.apply_translation((p + q) / 2.0)
        meshes.append(cyl)

    result = trimesh.util.concatenate(meshes)
    out_file = out_dir / f"max_volume_421_to_4d_{view['name']}.stl"
    result.export(out_file)
    print(f"Wrote {out_file}")
    print(f"  {len(result.faces)} triangles, {len(result.vertices)} mesh vertices")
