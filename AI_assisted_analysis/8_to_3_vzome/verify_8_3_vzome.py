#!/usr/bin/env python3
"""
STRICT SHAPE-RATIO VERIFICATION of the vZome 8->3 max-volume zonotope.

Pipeline A: Exact h from dV/dh=0 (proved identical to Q(rho) formula to
            50 decimal places via mpmath), standard Euclidean generators,
            all 256 projected 8-cube vertices.
Pipeline B: Parse vZome file, decode via heptagonal physical embedding,
            all 256 vertices.
Compare z_range / xy_diameter.
"""
import math, numpy as np, itertools, xml.etree.ElementTree as ET

print("=" * 70)
print("SHAPE RATIO VERIFICATION")
print("=" * 70)

# =====================================================================
# PIPELINE A: Exact Euclidean
# h^2 = (2/21)(1 - rho + rho^2) -- proved exact to 50 digits
# =====================================================================
print("\n--- PIPELINE A: Standard Euclidean (exact h from Q(rho)) ---\n")

rho_hp = 2 * math.cos(math.pi / 7)
h = math.sqrt(2/21 * (1 - rho_hp + rho_hp**2))
r = 4 / math.sqrt(21)
b = math.sqrt(8/3 - 7*h**2)

print(f"  r = {r:.17e}")
print(f"  h = {h:.17e}")
print(f"  b = {b:.17e}")

gens = []
for k in range(7):
    phi = 2 * math.pi * k / 7
    gens.append(np.array([r*math.cos(phi), r*math.sin(phi), h]))
gens.append(np.array([0.0, 0.0, b]))

verts_A = []
for signs in itertools.product([-1, 1], repeat=8):
    v = sum(s*g for s, g in zip(signs, gens))
    verts_A.append(v)
verts_A = np.array(verts_A)

z_range_A = np.max(verts_A[:,2]) - np.min(verts_A[:,2])
xy_radius_A = np.max(np.sqrt(verts_A[:,0]**2 + verts_A[:,1]**2))
xy_diam_A = 2 * xy_radius_A
ratio_A = z_range_A / xy_diam_A

print(f"  z_range     = {z_range_A:.17e}")
print(f"  xy_radius   = {xy_radius_A:.17e}")
print(f"  xy_diameter = {xy_diam_A:.17e}")
print(f"  RATIO = {ratio_A:.17e}")

# =====================================================================
# PIPELINE B: Parse vZome file
# =====================================================================
print("\n--- PIPELINE B: vZome file (heptagonal embedding) ---\n")

rho_val = 2 * math.cos(math.pi / 7)
sigma_val = rho_val**2 - 1
sigmaX2 = 2 * sigma_val
skewFactor = math.sin(3.0 / 7.0 * math.pi)

def decode(pt_str):
    v = list(map(int, pt_str.strip().split()))
    x_raw = v[0] + v[1]*rho_val + v[2]*sigma_val
    y_raw = v[3] + v[4]*rho_val + v[5]*sigma_val
    z_raw = v[6] + v[7]*rho_val + v[8]*sigma_val
    return np.array([x_raw + y_raw/sigmaX2, y_raw*skewFactor, z_raw])

vzome_path = (r"C:\Users\nanma\Documents\GitHub\hypercube-projections"
              r"\vZome\max-volume-8cube-to-3d.vZome")
tree = ET.parse(vzome_path)
root = tree.getroot()
for child in root:
    if "EditHistory" in child.tag:
        eh = child
        break

pts_B = []
for elem in eh:
    tag = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
    if tag == "ShowPoint":
        pts_B.append(decode(elem.get("point")))
pts_B = np.array(pts_B)

z_range_B = np.max(pts_B[:,2]) - np.min(pts_B[:,2])
xy_radius_B = np.max(np.sqrt(pts_B[:,0]**2 + pts_B[:,1]**2))
xy_diam_B = 2 * xy_radius_B
ratio_B = z_range_B / xy_diam_B

print(f"  Vertices    = {len(pts_B)}")
print(f"  z_range     = {z_range_B:.17e}")
print(f"  xy_radius   = {xy_radius_B:.17e}")
print(f"  xy_diameter = {xy_diam_B:.17e}")
print(f"  RATIO = {ratio_B:.17e}")

# =====================================================================
# COMPARISON
# =====================================================================
print("\n" + "=" * 70)
print("COMPARISON")
print("=" * 70)

abs_diff = abs(ratio_A - ratio_B)
rel_diff = abs_diff / ratio_A

print(f"\n  Pipeline A ratio: {ratio_A:.17e}")
print(f"  Pipeline B ratio: {ratio_B:.17e}")
print(f"  Absolute diff:    {abs_diff:.4e}")
print(f"  Relative diff:    {rel_diff:.4e}")

if rel_diff < 1e-13:
    verdict = "EXACT MATCH (within float64 precision)"
elif rel_diff < 1e-6:
    verdict = f"CLOSE but NOT exact (rel diff = {rel_diff:.2e})"
else:
    verdict = f"MISMATCH (rel diff = {rel_diff:.2e})"
print(f"\n  VERDICT: {verdict}")

scale_z = z_range_B / z_range_A
scale_xy = xy_diam_B / xy_diam_A
print(f"\n  Overall scale (z):  {scale_z:.17e}")
print(f"  Overall scale (xy): {scale_xy:.17e}")
print(f"  Scale ratio:        {scale_z/scale_xy:.17e} (should be 1.0)")
