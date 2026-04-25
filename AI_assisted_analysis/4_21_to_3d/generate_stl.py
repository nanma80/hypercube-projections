"""Generate STL of the max-volume 4_21->3D hull as balls (vertices) + cylinders (edges).

Outputs:
  max_volume_421_to_3d_balls_sticks.stl     (balls + cylinders)
  max_volume_421_to_3d_solid.stl            (the solid convex hull surface, for reference)
"""
import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from pathlib import Path

s2  = np.sqrt(2.0)
s23 = np.sqrt(2.0/3.0)   # = sqrt(6)/3
s21 = np.sqrt(21.0)

# 18 closed-form vertices in canonical orientation
V = []
# Hexagon at z=0, radius sqrt(2)
for i in range(6):
    a = i * np.pi/3
    V.append([s2*np.cos(a), s2*np.sin(a), 0.0])
# Top apex tri (z = +5/sqrt(21), angles 30,150,270)
for i in range(3):
    a = i*2*np.pi/3 + np.pi/6
    V.append([s23*np.cos(a), s23*np.sin(a), +5.0/s21])
# Upper mid (z = +4/sqrt(21), angles 90,210,330)
for i in range(3):
    a = i*2*np.pi/3 + np.pi/2
    V.append([s23*np.cos(a), s23*np.sin(a), +4.0/s21])
# Lower mid (z = -4/sqrt(21), angles 30,150,270)
for i in range(3):
    a = i*2*np.pi/3 + np.pi/6
    V.append([s23*np.cos(a), s23*np.sin(a), -4.0/s21])
# Bottom apex (z = -5/sqrt(21), angles 90,210,330)
for i in range(3):
    a = i*2*np.pi/3 + np.pi/2
    V.append([s23*np.cos(a), s23*np.sin(a), -5.0/s21])
V = np.array(V)
print(f"Vertices: {len(V)}")

# Get hull edges (the 48 unique edges)
hull = ConvexHull(V)
edge_set = set()
for simp in hull.simplices:
    for a, b in [(simp[0], simp[1]), (simp[1], simp[2]), (simp[0], simp[2])]:
        edge_set.add(tuple(sorted([a, b])))
edges = sorted(edge_set)
print(f"Hull volume: {hull.volume:.6f}  (target 8 sqrt(7)/3 = {8*np.sqrt(7)/3:.6f})")
print(f"Edges: {len(edges)}")

# Build balls + sticks
ball_radius = 0.08
stick_radius = 0.025

meshes = []
# Vertex balls
for v in V:
    s = trimesh.creation.icosphere(subdivisions=2, radius=ball_radius)
    s.apply_translation(v)
    meshes.append(s)
# Edge cylinders
for a, b in edges:
    p, q = V[a], V[b]
    seg = q - p
    L = np.linalg.norm(seg)
    cyl = trimesh.creation.cylinder(radius=stick_radius, height=L, sections=16)
    # cylinder is along z-axis by default at origin; align with seg.
    # Compute rotation from +z to seg direction.
    z = np.array([0., 0., 1.])
    direction = seg / L
    if np.allclose(direction, z):
        R = np.eye(4)
    elif np.allclose(direction, -z):
        R = trimesh.transformations.rotation_matrix(np.pi, [1, 0, 0])
    else:
        axis = np.cross(z, direction)
        axis /= np.linalg.norm(axis)
        angle = np.arccos(np.clip(np.dot(z, direction), -1, 1))
        R = trimesh.transformations.rotation_matrix(angle, axis)
    cyl.apply_transform(R)
    cyl.apply_translation((p + q) / 2.0)
    meshes.append(cyl)

balls_sticks = trimesh.util.concatenate(meshes)
out_dir = Path(__file__).parent
out1 = out_dir / "max_volume_421_to_3d_balls_sticks.stl"
balls_sticks.export(out1)
print(f"Wrote {out1}  ({len(balls_sticks.faces)} triangles)")

# Also export the solid convex hull surface
solid = trimesh.Trimesh(vertices=V, faces=hull.simplices, process=True)
out2 = out_dir / "max_volume_421_to_3d_solid.stl"
solid.export(out2)
print(f"Wrote {out2}  ({len(solid.faces)} triangles)")
