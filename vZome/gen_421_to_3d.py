"""
Generate vZome file for the 4_21 polytope (E8 root system, 240 roots)
projected to 3D using the SAME 8->3 projection that gives the
max-volume D_7 x Z_2 zonotope shadow of the 8-cube.

Same projection pi: R^8 -> R^3 (same 8 generator vectors); different input
vertices: 240 E8 roots instead of 256 cube vertices.

Result: a centrally-symmetric 3D point cloud with D_7d symmetry (order 28).
All coordinates exact in vZome's heptagonal field Q(rho).
"""
import numpy as np
import itertools

# 8 generator vectors g_1..g_8 in vZome's heptagonal-field encoding
# (1, rho, sigma) per coord, 9 ints per generator.
# These are 2x the natural generators (already scaled to clear denominators)
# from gen_vzome2.py / final_generators.py.
generators_vzome = [
    [-4, -2, -2,  2,  0,  2,  1, 2, 1],
    [ 2, -2, -4,  0,  0, -2,  1, 2, 1],
    [-2, -2,  2, -4, -2, -2,  1, 2, 1],
    [ 2,  4,  0,  2, -2, -4,  1, 2, 1],
    [ 0,  2,  4, -2, -2,  2,  1, 2, 1],
    [ 2,  0,  2,  2,  4,  0,  1, 2, 1],
    [ 0,  0, -2,  0,  2,  4,  1, 2, 1],
    [ 0,  0,  0,  0,  0,  0,  5, 4, 1],
]

rho_val = 2*np.cos(np.pi/7)
sigma_val = rho_val**2 - 1
sigmaX2 = 2 * sigma_val
skewFactor = np.sin(3.0/7.0 * np.pi)

def vzome_to_phys(v):
    x = v[0] + v[1]*rho_val + v[2]*sigma_val
    y = v[3] + v[4]*rho_val + v[5]*sigma_val
    z = v[6] + v[7]*rho_val + v[8]*sigma_val
    return np.array([x + y/sigmaX2, y*skewFactor, z])

def vadd(a, b):  return [a[i]+b[i] for i in range(9)]
def vsub(a, b):  return [a[i]-b[i] for i in range(9)]
def vneg(a):     return [-a[i] for i in range(9)]
def vscale(a,s): return [a[i]*s for i in range(9)]

# Build the 240 E8 roots and project each to a vZome 9-tuple.
# To keep everything integer, project 2*root:
#   - integer root r = e_i +/- e_j  ==>  2r projects to 2*(g_i +/- g_j)
#   - half-integer root r = (1/2)(s_1,...,s_8)  ==>  2r projects to sum_i s_i g_i

projected = []  # list of (root_label, vzome_9tuple)

# 112 integer roots: signs (e_i, e_j) with i<j and 4 sign combos
int_count = 0
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product([+1,-1], repeat=2):
        v = [0]*9
        v = vadd(v, vscale(generators_vzome[i], 2*si))
        v = vadd(v, vscale(generators_vzome[j], 2*sj))
        projected.append((f"int_{i}{j}_{si:+d}{sj:+d}", v))
        int_count += 1
assert int_count == 112

# 128 half-integer roots: (s_1,...,s_8) with prod = +1
hi_count = 0
for signs in itertools.product([+1,-1], repeat=8):
    if np.prod(signs) != 1: continue
    v = [0]*9
    for k in range(8):
        v = vadd(v, vscale(generators_vzome[k], signs[k]))
    projected.append((f"hi_{''.join('+' if s>0 else '-' for s in signs)}", v))
    hi_count += 1
assert hi_count == 128
assert len(projected) == 240

# Deduplicate identical projected points
unique = {}
for label, v in projected:
    key = tuple(v)
    unique.setdefault(key, []).append(label)

print(f"Total roots:           {len(projected)}")
print(f"Distinct image points: {len(unique)}")
mult = sorted({len(v) for v in unique.values()})
for m in mult:
    n = sum(1 for v in unique.values() if len(v)==m)
    print(f"  multiplicity {m}: {n} image points  (covers {m*n} roots)")

# Physical positions (for sanity)
phys = np.array([vzome_to_phys(list(k)) for k in unique.keys()])
print(f"\nBounding box (physical, units of half-generators):")
print(f"  x: [{phys[:,0].min():+.3f}, {phys[:,0].max():+.3f}]")
print(f"  y: [{phys[:,1].min():+.3f}, {phys[:,1].max():+.3f}]")
print(f"  z: [{phys[:,2].min():+.3f}, {phys[:,2].max():+.3f}]")

# ============================================================
# Verify D_7d symmetry using induced 3D matrices from W(E8) elements
# (the C_7 cyclic permutation, a transposition reflection, and -I).
# ============================================================
G_phys_gens = np.array([vzome_to_phys(g) for g in generators_vzome])
def induced_3d(perm):
    """3D matrix induced by a permutation of {g_0,..,g_7} via least-squares."""
    return G_phys_gens[perm].T @ np.linalg.pinv(G_phys_gens.T)

theta = 2*np.pi/7
C7_3d = np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
T_refl = induced_3d([1,0,6,5,4,3,2,7])  # reflection from transposition (0,1) in S_7

def max_isometry_residual(M):
    transformed = phys @ M.T
    return max(np.min(np.linalg.norm(phys - tp, axis=1)) for tp in transformed)

print(f"\nC_7 rotation residual:           {max_isometry_residual(C7_3d):.2e}")
print(f"D_7 reflection residual:         {max_isometry_residual(T_refl):.2e}")
print(f"Central inversion -I residual:   {max_isometry_residual(-np.eye(3)):.2e}")
# Enumerate generated group
gens_3d = [C7_3d, T_refl, -np.eye(3)]
group=[np.eye(3)]
for _ in range(6):
    new=[]
    for g in group:
        for h in gens_3d:
            c = g @ h
            if not any(np.allclose(c,x,atol=1e-8) for x in group+new):
                new.append(c)
    if not new: break
    group += new
preserving = [g for g in group if max_isometry_residual(g) < 1e-6]
print(f"Symmetry group order: {len(preserving)} (expected D_7d = 28)")

# ============================================================
# Write vZome file
# ============================================================
def pt_str(v):
    return " ".join(str(int(x)) for x in v)

cmds = []
for key in unique.keys():
    cmds.append(f'    <ShowPoint point="{pt_str(key)}"/>')

# Edges: connect each E8 root to its closest neighbors using E8's
# minimum-distance edge structure. Two roots r1, r2 are joined by an edge
# of E8 (in the 8-dim sense, i.e. <r1,r2> = 1) iff r1 - r2 is itself a root.
# We use this to draw edges in the projection.
# Build set of root vectors in R^8 (multiplied by 2 to keep ints):
roots_8d = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product([+1,-1], repeat=2):
        r = [0]*8; r[i] = 2*si; r[j] = 2*sj  # 2*root for integer type
        roots_8d.append(tuple(r))
for signs in itertools.product([+1,-1], repeat=8):
    if np.prod(signs) != 1: continue
    roots_8d.append(tuple(signs))  # 2*(half-integer root)
assert len(roots_8d) == 240

# Map root_8d (tuple) -> projected vzome key (tuple)
root_to_proj = {}
idx = 0
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product([+1,-1], repeat=2):
        r = [0]*8; r[i] = 2*si; r[j] = 2*sj
        key = tuple(projected[idx][1]); root_to_proj[tuple(r)] = key; idx += 1
for signs in itertools.product([+1,-1], repeat=8):
    if np.prod(signs) != 1: continue
    key = tuple(projected[idx][1]); root_to_proj[tuple(signs)] = key; idx += 1

root_set = set(roots_8d)
edges_set = set()
for r1 in roots_8d:
    for r2 in roots_8d:
        if r1 >= r2: continue
        diff = tuple(r1[k]-r2[k] for k in range(8))
        # diff (or its half) must be a root in 4*-scaled or 2*-scaled form
        # Since r1, r2 are 2*roots, diff is 2*(root1-root2). For E8, root1-root2
        # is a root iff <root1,root2> = 1.
        # <root1,root2> = (1/4) <r1,r2> in 2x scaling.
        ip = sum(r1[k]*r2[k] for k in range(8))
        # original inner product = ip/4. We want it = 1 for adjacency.
        if ip == 4:
            p1, p2 = root_to_proj[r1], root_to_proj[r2]
            if p1 != p2:
                edges_set.add((min(p1,p2), max(p1,p2)))

print(f"\nE8 edges (root pairs with <r1,r2>=1): {sum(1 for r1 in roots_8d for r2 in roots_8d if r1<r2 and sum(r1[k]*r2[k] for k in range(8))==4)}")
print(f"Distinct projected edges: {len(edges_set)}")

for p1, p2 in sorted(edges_set):
    cmds.append(f'    <JoinPointPair start="{pt_str(list(p1))}" end="{pt_str(list(p2))}"/>')

xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<vzome:vZome xmlns:vzome="http://xml.vzome.com/vZome/4.0.0/"
             buildNumber="111" field="heptagon" version="7.1">
  <EditHistory editNumber="{len(cmds)}" lastStickyEdit="-1">
{chr(10).join(cmds)}
  </EditHistory>
  <notes/>
  <sceneModel ambientLight="41,41,41" background="175,200,220">
    <directionalLight color="235,235,228" x="1.0" y="-1.0" z="-1.0"/>
    <directionalLight color="228,228,235" x="-1.0" y="0.0" z="0.0"/>
    <directionalLight color="30,30,30" x="0.0" y="0.0" z="-1.0"/>
  </sceneModel>
  <Viewing>
    <ViewModel distance="100.0" far="200.0" near="0.25"
               parallel="false" stereoAngle="0.0" width="45.0">
      <LookAtPoint x="0.0" y="0.0" z="0.0"/>
      <UpDirection x="0.0" y="1.0" z="0.0"/>
      <LookDirection x="0.0" y="0.0" z="-1.0"/>
    </ViewModel>
  </Viewing>
  <SymmetrySystem name="heptagonal antiprism corrected" renderingStyle="heptagonal antiprism">
    <Direction color="217,18,24" name="red" orbit="red"/>
    <Direction color="0,153,63" name="green" orbit="[[-2,0,1,1],[0,0,0,1]]"/>
    <Direction color="0,142,194" name="blue" orbit="[[0,0,0,1],[0,0,0,1]]"/>
  </SymmetrySystem>
  <OtherSymmetries/>
  <Tools/>
</vzome:vZome>
'''

outpath = r"C:\Users\nanma\Documents\GitHub\hypercube-projections\vZome\heptagonal-421-to-3d.vZome"
with open(outpath, 'w', encoding='utf-8') as f:
    f.write(xml)
print(f"\nWrote: {outpath}")
print(f"  {len(unique)} distinct vertices, {len(edges_set)} edges, {len(cmds)} edits")
