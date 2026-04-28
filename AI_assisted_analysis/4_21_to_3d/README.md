# Maximum-Volume 3D Shadow of the 4_21 Polytope

## Setup

The **4_21 polytope** has the 240 E₈ roots as its vertices in ℝ⁸:
- **112 integer roots:** all signed permutations of (±1, ±1, 0, 0, 0, 0, 0, 0)
- **128 half-integer roots:** (±½, ±½, ±½, ±½, ±½, ±½, ±½, ±½) with an even
  number of minus signs

Every root has Euclidean norm √2.

We seek the 3-dimensional subspace W ⊂ ℝ⁸ that maximizes
**Vol₃(proj_W(4_21))** — i.e. the convex-hull volume of the 240 projected
points.

The convex-hull volume cannot be written as a sum of |det| over k-subsets
(4_21 is **not** a zonotope), so we use `scipy.spatial.ConvexHull` and
basin-hopping over ℝ³ˣ⁸ matrices.

## Result

**Max volume = 8√7 / 3 ≈ 7.0553368295**, computed with the
basis-orthonormalization convention (each E₈ root has length √2 in ℝ⁸,
projection basis is orthonormal). Reproducible with seeds 7 and 42;
30 basin-hopping trials all converged to this value (a few stuck at
the slightly lower local maximum 7.04630881). The closed form
8√7/3 was identified analytically (see "Closed-form solution" below)
and verified to ∼15 digits.

| Projection                           | Vol      | Ratio  |
|--------------------------------------|----------|--------|
| **Max-volume (this work)**           | 7.05534  | 1.000  |
| 7-fold heptagonal (D₇d, vZome model) | 6.26678  | 0.888  |

So the heptagonal symmetric projection achieves about **89%** of the maximum.

### Hull combinatorics

| | Max-area (2D) | Max-volume (3D) | 7-fold (D₇d, 3D) | 4D max-volume | H₄ to 4D (600-cell) |
|---|---|---|---|---|---|
| Vertices | **6** | **18** (3 orbits of 6) | 240 | 48 | 120 |
| Edges    | **6** | **48** (4 orbits of 12) | (huge) | 144 | 720 |
| Facets   | — | **32** (all triangles) | (huge) | n/a | n/a |
| Symmetry | **D₆, order 12** | **D₃d, order 12** (with −I) | D₇d, order 28 | unknown | H₄, order 14400 |
| Equilateral hull verts? | Yes | No (3 norm classes) | Yes | No (2 norms) | Yes |
| Volume (orthonormal frame) | 3√3 ≈ 5.196 | 8√7/3 ≈ 7.055 | 6.267 | 8.880 | 8.090 |

**Note:** The H₄ projection produces a **600-cell** with **120** vertices (not 600 —
the "600" refers to tetrahedral cells). The original `e8.py` comments stated
"600 vertices", which was an error corrected in the 4D re-verification. Volumes
in the old `e8.py` (129.44, 142.08) used vertices scaled by 2×, so they are
16× larger than the unit-root values shown here.

**The shape is not a named polyhedron.** It has 4 distinct edge lengths
(≈ 0.8452, 1.1952, 1.3628, 1.4142), so it is not one of the 8 convex
deltahedra; its (V, E, F) = (18, 48, 32) and degree sequence
(12 hexavalent + 6 tetravalent vertices) do not match any Johnson,
Archimedean, or Catalan solid. It has trigonal antiprismatic symmetry
D₃d (order 12) — same Schoenflies family as the heptagonal model
(D₇d), but with a 3-fold instead of 7-fold principal axis.

The 3 vertex orbits sit at norms √2 (those 6 are E₈ roots that survive
the projection at full length), 1.3628, and 1.1952 (vertices that get
shortened by projection).

The previous 4D study in `e8.py` already revealed that the H₄-symmetric
projection (whose hull is the 600-cell) is **not** the 4D max-volume
projection — the 4D max gives 48 vertices and is non-equilateral. The
3D situation is analogous: the max-volume hull is non-equilateral and
has only 18 vertices.

## Striking structural finding: 4 axes coalesce

The optimal 3×8 projection matrix Q (orthonormal columns) consistently
has **four identical rows**, e.g. (different runs may permute which
4 of the 8 indices are involved):

    row a = row b = row c = row d ≈ (0.1565, -0.0989, 0.4547)

In other words: **4 of the 8 ℝ⁸ axes project to a single direction in ℝ³.**

This means the kernel of the optimal 8→3 map contains the 3-dimensional
subspace
    K₀ = span{ e_a − e_b, e_b − e_c, e_c − e_d }
inside the 4-dimensional coordinate subspace {axes a, b, c, d}. The full
kernel is 5-dimensional (since the projection has rank 3); K₀ accounts
for 3 of those 5 dimensions, the other 2 are "generic".

## Factorization 8 → 5 → 3 (answers user question)

**Yes**: because K₀ ⊂ ker(π), the optimal 3D projection π : ℝ⁸ → ℝ³
factors through ℝ⁸/K₀ ≅ ℝ⁵:

    ℝ⁸  ──A──►  ℝ⁵  ──B──►  ℝ³

where A is the **explicit "coalesce 4 axes"** linear map (R⁵-orthonormal
columns):

    column 0 of A = (½, ½, ½, ½, 0, 0, 0, 0)ᵀ      ← coalesced 4 axes
    column 1..4 of A = e₅, e₆, e₇, e₈              ← kept axes

and B is some 3D subspace of ℝ⁵.

**Empirical verification** (`max_shadow_421_to_3d.py` + the inline
factorization check): performing this 8→5 step *first*, then maximizing
the 5→3 convex-hull volume over the resulting 83-point set (240 4_21
vertices collapse to **83 distinct points** in ℝ⁵), recovers the same
maximum **7.0553368295** to 10 decimal places. So no volume is lost by
restricting to the 8→5→3 factorization; the 5D intermediate point set
is sufficient.

A few corollaries / caveats:

- This factorization is **not** automatic for arbitrary projections.
  Trivially every linear ℝ⁸→ℝ³ map factors through any intermediate ℝ⁵,
  but only because we *choose* the intermediate ℝ⁵ to contain the kernel.
  The interesting fact is that **at the optimum the kernel happens to
  contain the highly symmetric subspace K₀** (a 3-dim subspace where 4
  coordinate axes all become equal). That is a structural property of
  the 4_21 maximum, not a generic feature of optimization.
- Equivalently: the optimum sits on the symmetric stratum where the
  symmetric group S₄ acting by permuting the 4 selected axes is in the
  stabilizer of the projection (because all 4 axes are sent to the same
  3D direction, any permutation of them leaves Q unchanged).
- The 5D image of 4_21 under A is a 240→83 point set whose coordinates
  all lie in {−1, −½, 0, +½, +1}. Identifying this 5D polytope by name
  (if any) is left as a follow-up.

The lower (spurious) local max at 7.04631 corresponds to a different
coalescence pattern (likely 3+2 or 5+1 instead of 4+1+1+1+1) and is
not the global maximum.

## Files

- `../../max_shadow_421_to_3d.py` — main script. Builds the 240 E₈ roots,
  computes the heptagonal baseline, and runs basin-hopping over ℝ³ˣ⁸
  matrices to find the max-volume 3D shadow.
- `../../viz_421_to_3d.py` — finds the optimum, canonicalizes orientation
  (C₃ axis aligned with +z, one C₂ axis with +x), saves generators, and
  produces a 3D matplotlib visualization.
- `max_volume_421_to_3d_generators.txt` — the 8 canonicalized generators
  g₀…g₇ in ℝ³, plus the 18 hull-vertex coordinates layered by z.
- `max_volume_421_to_3d_generators.png` — single-view 3D plot.
- `max_volume_421_to_3d_views.png` — three-view plot (perspective, top-down,
  side) showing the 5-layer structure (3 + 3 + 6 + 3 + 3 vertices along
  the C₃ axis).
- `../../e8.py` — earlier 4_21 → 4D study. Found 4D max ≈ 142.08 with
  48 hull vertices (vs 129.44 with 600 vertices for the H₄ projection).
- `../../vZome/heptagonal-421-to-3d.vZome` — the 7-fold (D₇d) model
  in vZome, used as the 3D baseline above.
- `../../vZome/gen_421_to_3d.py` — generator for that vZome model.
- `generate_stl.py` — builds STL meshes of the max-volume 18-vertex hull
  (vertices as balls, edges as cylinders) and the bare convex-hull
  surface, using `trimesh`.
- `max_volume_421_to_3d_balls_sticks.stl` — ball-and-stick STL of the
  hull only (18 vertex spheres + 48 edge cylinders, all uniform size).
- `max_volume_421_to_3d_solid.stl` — the bare 32-triangle convex-hull
  surface as STL.
- `generate_stl_all_struts.py` — builds the richer "all 4_21 struts
  projected into the hull interior" STL (49 balls + 432 cylinders),
  with sizes scaled by **multiplicity** (see next section).
- `max_volume_421_to_3d_all_struts.stl` — STL of the full projected
  4_21 strut graph: every distinct projected vertex and every distinct
  projected edge, including those interior to the hull.

### Visualizing all projected 4_21 struts (`*_all_struts.stl`)

The 4_21 polytope has **240 vertices and 6720 edges** (one for each
pair of E₈ roots with dot product 1). Under the max-volume projection,
many of these collapse:

- **240 → 49 distinct 3D points**: 18 of them are the hull vertices,
  the other 31 lie *inside* the hull.
- **6720 edges → 432 distinct nonzero 3D segments** (672 edges
  collapse to zero length).

Each 3D point and each 3D segment has a well-defined **multiplicity**
= the number of 8D objects that map to it. We use the multiplicity to
set a small number of discrete sizes, so the visualization shows
which points / struts are "richer".

**Vertex (ball) multiplicities — three values, sizes 1.0× / 1.5× / 2.0×
of base radius:**

| 8D root multiplicity | # of 3D points | Ball size | Where               |
|----------------------|----------------|-----------|---------------------|
| 1                    | 24             | 1.0×      | all 18 hull verts + 6 interior |
| 8                    | 24             | 1.5×      | interior            |
| 24                   | 1              | 2.0×      | the **origin**      |

(Sum check: 24·1 + 24·8 + 1·24 = 240 ✓.)

**Strut (cylinder) multiplicities — four values, sizes 1.0× / 1.5× /
2.0× / 2.5× of base radius:**

| 8D edge multiplicity | # of 3D struts | Strut size | 8D edges contributed |
|----------------------|----------------|------------|----------------------|
| 1                    | 96             | 1.0×       | 96                   |
| 8                    | 216            | 1.5×       | 1728                 |
| 32                   | 96             | 2.0×       | 3072                 |
| 48                   | 24             | 2.5×       | 1152                 |

(Sum check: 96·1 + 216·8 + 96·32 + 24·48 = 6048 = 6720 − 672 ✓.)

**Striking pattern.** Every strut whose length matches one of the four
hull-edge lengths (√(5/7), √(10/7), √(13/7), √2) decomposes as
**24 mult-1 + 18 mult-8** struts. The five interior-only lengths
(0.327, 0.655, 0.824, 0.926, 0.982) have multiplicities only in
{8, 32, 48} — no strut interior to the hull is hit by a single E₈
edge. So mult-1 struts are exactly those whose endpoints are both
hull-orbit (mult-1) vertices. The biggest "highway" struts (mult 48)
all pass through the origin or near-origin region.

The constants `base_ball_r` and `base_stick_r` at the top of
`generate_stl_all_struts.py` control overall thickness; the per-bin
multipliers `(1.0, 1.5, 2.0, 2.5)` control the relative emphasis.

### Why the max-volume shape is **not** constructible exactly in vZome

The 18 hull-vertex coordinates live in **Q(√2, √3, √7)** (degree 8 over Q):
- x ∈ {0, ±1/√2, ±√2}            → needs √2
- y ∈ {0, ±1/√6, ±√6/3, ±√6/2}    → needs √6 = √2·√3
- z ∈ {0, ±4/√21, ±5/√21}         → needs √21 = √3·√7

The killer is **√7**. None of vZome's standard fields contains √7:

| Field        | Generator        | Contains √7? |
|--------------|------------------|--------------|
| Golden       | √5 (φ)           | no           |
| Root-2       | √2               | no           |
| Root-3       | √3               | no           |
| Heptagonal   | 2·cos(π/7) (degree 3) | **no**  |
| SnubDodec, Plastic, SqrtPhi, …                  | no           |

The heptagonal field — which one might naively expect to "contain √7" —
in fact does not. It is the maximal real subfield of Q(ζ₇), and the only
quadratic subfield of Q(ζ₇) is Q(√(−7)) (because 7 ≡ 3 mod 4), which is
purely imaginary.

Crucially, the field requirement is **isometry-invariant**: rotating the
shape into a different orientation cannot eliminate any of √2, √3, √7
from the coordinate set. (All squared inter-vertex distances are
rationals with denominator 7, and the rank-3 rational Gram matrix
forces realization in Q(√d₁, √d₂, √d₃) with Cholesky pivots 2, 6, 21.)
The ⁄7 flavour is baked into the maximum because

> V_max² = 448/9 = 64·7/9,  hence V_max = 8√7/3.

So a faithful vZome model would require a Q(√2, √3, √7) field that does
not currently exist. Practical fall-backs (e.g. anisotropic z-stretch
to put z in Q while keeping xy in Q(√2,√3); or numerical placement in
the golden field) preserve only combinatorics or only approximate
geometry. For exact 3D visualization we therefore use STL instead
(see files above).

### Construction recipe (for reproducing the shape)

The 8 canonicalized generators (rows of the 8×3 orthonormal matrix Q,
recorded numerically in `max_volume_421_to_3d_generators.txt`):

      g_0  ≈ (+1/√2,  0,       +0.16366)        (~half-norm 0.7258)
      g_1  ≈ (−1/√2,  0,       +0.16366)        (~half-norm 0.7258)
      g_2  ≈ ( 0,    −1/√6,    +0.38188)        (~half-norm 0.5590)
      g_3  ≈ ( 0,    +1/√6,    −0.70921)        (~half-norm 0.8183)
      g_4=g_5=g_6=g_7
           ≈ ( 0,    +1/√6,    +0.27277)        (~half-norm 0.4910)

(All but g₀, g₁ have x = 0, so they lie in the yz-plane. The four
"coalesced" axes 4–7 share a single direction.)

The 240 vertices in ℝ³ are then exactly:

      ±g_i ± g_j               for 0 ≤ i < j ≤ 7         (112 vertices)
      (1/2) Σ_{k=0..7} s_k g_k for s ∈ {±1}⁸ with even Π (128 vertices)

Hull vertices line up in 5 layers along the C₃ (=z) axis:

| z          | count | vertex norm | ring shape           |
|------------|-------|-------------|-----------------------|
| +1.0911    |   3   | 1.3628      | small triangle (top)  |
| +0.8729    |   3   | 1.1952      | small triangle, twisted 60° |
|  0         |   6   | √2 = 1.4142 | regular hexagon       |
| −0.8729    |   3   | 1.1952      | small triangle        |
| −1.0911    |   3   | 1.3628      | small triangle, twisted 60° |

The numerical constants 1/√2, 1/√6, √2, etc. that appear in the
canonicalized generators turn out to extend to a fully algebraic
optimum — see next section.

### Closed-form solution

Symbolically, the canonicalized generators are

      g_0  = (+1/√2,  0,       +3/(4√21))
      g_1  = (−1/√2,  0,       +3/(4√21))
      g_2  = ( 0,    −1/√6,    +7/(4√21))
      g_3  = ( 0,    +1/√6,   −13/(4√21))
      g_4=g_5=g_6=g_7
           = ( 0,    +1/√6,    +5/(4√21))

i.e. all four z-components share the denominator **4√21**, with
integer numerators **(3, 3, 7, −13, 5, 5, 5, 5)** that sum to zero
(centering condition), satisfy the y/z-orthogonality `−7 + (−13) + 4·5 = 0`,
and have squared sum `2·9 + 49 + 169 + 4·25 = 336 = 16·21` so that
each column of Q is a unit vector.

**Maximum volume:** plugging these exact generators into the 18-vertex
hull formula yields

> **V_max = 8√7 / 3**

(verify: V² = 448/9, matched to 15 digits by both basin-hopping and
direct evaluation of `ConvexHull` on the symbolic vertices).

**Two key parameters.** With `h = z_top = 4d` and `k = z_mid = a − c`,
the 18 hull vertex positions depend only on (h, k), and the orthonormality
of Q forces them onto the constraint ellipse

> **h² − hk + k² ≤ 1**.

The maximum lies on the boundary, at **h = 5/√21**, **k = 4/√21**
(equivalently the simple ratio **h : k = 5 : 4**), where the
discriminant of the quadratic for c collapses to zero and pins c
uniquely to **−13/(4√21)**.

The 18-vertex polyhedron in those coordinates:
- equatorial regular hexagon at z = 0, radius √2;
- four C₃-symmetric triangles at z = ±h = ±5/√21 (radius √(2/3),
  vertex angles 30°, 150°, 270°) and z = ±k = ±4/√21 (same radius,
  angles 90°, 210°, 330°).

### Closed-form vertex list (all 18)

In the canonical orientation (C₃ axis along +z, one C₂ axis along +x):

**Hexagon, z = 0** (vertex norm √2):

| label | x       | y      | z |
|-------|---------|--------|---|
| H₀    | +√2     | 0      | 0 |
| H₁    | +1/√2   | +√6/2  | 0 |
| H₂    | −1/√2   | +√6/2  | 0 |
| H₃    | −√2     | 0      | 0 |
| H₄    | −1/√2   | −√6/2  | 0 |
| H₅    | +1/√2   | −√6/2  | 0 |

**Top apex triangle, z = +5/√21** (vertex norm √(13/7)):

| label | x       | y       | z       |
|-------|---------|---------|---------|
| T₀    | +1/√2   | +1/√6   | +5/√21  |
| T₁    | −1/√2   | +1/√6   | +5/√21  |
| T₂    | 0       | −√6/3   | +5/√21  |

**Upper mid triangle, z = +4/√21** (vertex norm √(10/7)):

| label | x       | y       | z       |
|-------|---------|---------|---------|
| M₀    | 0       | +√6/3   | +4/√21  |
| M₁    | −1/√2   | −1/√6   | +4/√21  |
| M₂    | +1/√2   | −1/√6   | +4/√21  |

**Lower mid triangle, z = −4/√21** (vertex norm √(10/7)):

| label | x       | y       | z       |
|-------|---------|---------|---------|
| N₀    | +1/√2   | +1/√6   | −4/√21  |
| N₁    | −1/√2   | +1/√6   | −4/√21  |
| N₂    | 0       | −√6/3   | −4/√21  |

**Bottom apex triangle, z = −5/√21** (vertex norm √(13/7)):

| label | x       | y       | z       |
|-------|---------|---------|---------|
| B₀    | 0       | +√6/3   | −5/√21  |
| B₁    | −1/√2   | −1/√6   | −5/√21  |
| B₂    | +1/√2   | −1/√6   | −5/√21  |

Useful identities: `√(2/3) = √6/3`, `√6/2 = √(3/2)`. Central
inversion gives B = −T and N = −M; the mid layers are rotated 60°
relative to the apex layers. All x,y coordinates lie in
`{0, ±1/√2, ±1/√6, ±√6/3, ±√6/2, ±√2}`; all z values in
`{0, ±4/√21, ±5/√21}`.

### Geometric structure: a "squashed cuboctahedron" with six tents

A satisfying geometric description of the 18-vertex hull:

**Step 1 — drop the M, N orbit.** The 12 remaining vertices `{H, T, B}`
form a polytope that is **combinatorially a cuboctahedron**: 12 vertices,
24 edges, **8 triangles + 6 quadrilaterals**, with the same vertex-figure
incidence as the Archimedean cuboctahedron.

But it is **not metrically** a cuboctahedron:

| Quantity                         | Our `{H,T,B}` polytope    | True cuboctahedron |
|----------------------------------|---------------------------|--------------------|
| Hexagon radius                   | √2                        | √2                 |
| Apex-triangle xy-radius          | √(2/3)                    | √(2/3)             |
| Apex-triangle z-height           | **5/√21 ≈ 1.0911**        | 2/√3 ≈ 1.1547      |
| Vertex norms                     | {√2, √(13/7)}             | √2 (uniform)       |
| Hexagon and apex-tri edge        | √2                        | 1 (after rescaling)|
| Lateral H–T / H–B edge           | **√(13/7) ≠ √2**          | equal to other     |
| Six "square" faces               | **rectangles** √2 × √(13/7) | true squares     |

So `{H, T, B}` is a triangular gyrobicupola **squashed along its C₃ axis**
by the factor **5/(2√7) ≈ 0.945** relative to the regular Archimedean
cuboctahedron. Its quad faces are rectangles with aspect ratio
√(13/14) ≈ 0.964, not squares. (Direct check: the lateral edges of each
quad are perpendicular to both the hexagon edge and the apex-triangle
edge, confirming "rectangle", not "isosceles trapezoid".)

**Step 2 — add a tent over each rectangle.** Each of the 6 mid vertices
`M`/`N` sits **directly along the outward radial axis through the centre
of one rectangle face**, at radius bigger than the rectangle's centroid:

- `M₀ = (0, √6/3, 4/√21)` lies along the radial line of the rectangle
  (H₁, T₀, T₁, H₂) (whose centroid is (0, √6/3, 5/(2√21)), same x and y).
  Since 4/√21 > 5/(2√21), `M₀` is pushed *outward* of the rectangle.
- The other M, N vertices play the same role for the other five rectangle
  faces (3 covered by M's, 3 by N's, alternating around the equator).

So the max-volume hull is built by **erecting 6 four-sided pyramids
("tents") over the 6 rectangle faces of the squashed cuboctahedron**,
each tent's apex being one of the M / N vertices.

**Face/edge accounting through this lens:**

- 8 triangle faces of the cuboctahedron skeleton survive untouched
  (1 top cap + 1 bottom cap + 6 lateral H–H–T / H–H–B).
- Each of the 6 rectangles is replaced by **4 triangles** of the tent,
  giving 6 × 4 = 24 new triangle faces.
- Total: 8 + 24 = **32 triangles** ✓.
- Edges: 24 (cuboctahedron skeleton) + 4·6 = 24 (tent ridge-edges
  M–H, M–T or N–H, N–B) = **48** ✓.
- Vertices: 12 (skeleton) + 6 (tent apices) = **18** ✓.

This decomposition also makes the volume easy to organise as
`V = V_squashed_cubocta + 6·V_tent`.

### Closed-form edges and faces of the hull


> **18 vertices, 48 edges, 32 faces — all triangles.**

(Euler: 18 − 48 + 32 = 2 ✓.) All four distinct **edge lengths** are
square roots of rationals with denominator 7:

| Edge type                                            | Length     | Count |
|------------------------------------------------------|------------|-------|
| H–H (hexagon side) and T–T, B–B (top/bot tri sides)  | **√2 = √(14/7)** | 12 |
| H–T and H–B  (hexagon to top/bot apex tri)           | **√(13/7)**       | 12 |
| H–M and H–N  (hexagon to middle apex tri)            | **√(10/7)**       | 12 |
| T–M and B–N  (top to mid; bottom to mid)             | **√(5/7)**        | 12 |

The 32 triangular faces fall into **5 D₃d-orbits**:

| Face type (vertex labels)         | Orbit size | Description                                          |
|-----------------------------------|------------|------------------------------------------------------|
| (T,T,T) / (B,B,B)                 | 1+1 = 2    | the equilateral cap at the top and bottom            |
| (T,T,M) / (B,B,N)                 | 3+3 = 6    | between an apex-tri edge and a mid-tier vertex       |
| (M,T,T) ✗ already listed          |            |                                                      |
| (H,H,T) / (H,H,B)                 | 3+3 = 6    | upward/downward triangles, hexagon edge to apex tri  |
| (H,H,M) / (H,H,N)                 | 3+3 = 6    | hexagon edge to mid-tier vertex                      |
| (H,M,T) / (H,N,B)                 | 6+6 = 12   | skew triangles spanning all three latitudes          |

Total: 2 + 6 + 6 + 6 + 12 = 32. (No quadrilateral or pentagonal faces;
no merged coplanar triangles either — every triangulated facet is its
own flat face.)

**Vertex degrees:** the hexagonal vertices and the top/bot apex-tri
vertices all have degree 6 (each touches 6 triangles). The mid-tier
vertices `M, N` have degree 4. (Verified: 6·6 + 6·6 + 6·4 = 96 = 2·48.)

### Methodology: numerical-first, then conjecture-and-verify

The closed form was *not* derived from first principles. The flow
that produced it was:

1. **Numerical maximum.** Run basin-hopping over the full 24-parameter
   3×8 projection matrix (`max_shadow_421_to_3d.py`) → V ≈ 7.05533683.
2. **Observe structure.** Canonicalize the optimum
   (`viz_421_to_3d.py`) and inspect: x,y components look like
   1/√2, 1/√6; four input axes coalesce; the symmetry group is D₃d.
3. **Conjecture a parametric ansatz.** Freeze the x,y components at
   the recognized exact values, leaving only four z-components
   (a, b, c, d). Apply orthonormality of Q to cut two more parameters
   (`−b + c + 4d = 0` and the unit-norm constraint), reducing to
   **two free parameters** (h, k).
4. **Re-optimize within the ansatz.** Maximize on the feasibility
   ellipse `h² − h k + k² ≤ 1` and confirm the volume matches the
   numerical max — this confirms the ansatz captures the optimum.
5. **Spot rational signatures.** In the 2-parameter optimum,
   `k/h = 0.80000…`, `h² = 1.19047… ≈ 25/21`, `(h−k)/(h+k) = 0.11111…`
   → conjecture **h = 5/√21, k = 4/√21** (h : k = 5 : 4).
6. **Push through symbolically.** With those h, k the c-quadratic's
   discriminant is exactly zero, so c is pinned to −13/(4√21); the
   z-components turn out to be integers (3, 7, −13, 5)/(4√21).
7. **Recognize the volume.** V² ≈ 49.7777… = 448/9 → conjecture
   **V = 8√7/3**.
8. **Verify.** Plug the symbolic generators back through ConvexHull
   in 50-digit mpmath; the volume agrees with 8√7/3 to 15 digits and
   Q is exactly orthonormal.

Strictly speaking this proves only that 8√7/3 is the value at a
specific D₃d critical point, not that it is the global maximum over
all 3D subspaces. But it agrees with every basin-hopping run, so the
conjecture-and-verify chain is strong evidence. A full proof would
require ruling out other symmetry classes and hull combinatorial
types — see open questions below.

## Open questions

1. ~~What is the **symmetry group** of the max-volume 3D shadow?~~
   **Answer: D₃d (trigonal antiprismatic), order 12**, computed by
   brute-force enumeration of orthogonal matrices preserving the
   18-vertex set. The group decomposes as 1 identity + 2 C₃ rotations
   + 3 C₂ rotations + 6 improper (including −I). The shape itself is
   not a named polyhedron.
2. ~~Does the volume `7.05533683…` admit an algebraic closed form?~~
   **Answer: yes, V_max = 8√7/3.** The optimum has h = 5/√21,
   k = 4/√21, and integer-numerator z-components (3, 7, −13, 5)/(4√21).
3. Is the 5D intermediate (83-point set) a recognizable polytope?
4. ~~Does this 4-axis-coalescence pattern generalize? (E.g. for 4_21 → 4D,
   do some axes coalesce at the optimum? `e8.py` did not check.)~~
   **Answer: No.** Re-verification with `max_shadow_421_to_4d.py` confirms
   that the 4D max-volume projection has **no generator coalescence** — all
   8 generators are distinct (though they split into 2 norm classes of 4).
   The coalescence pattern is dimension-dependent:
   - 2D max: massive coalescence (6 of 8 axes coalesce in pairs)
   - 3D max: 4 axes coalesce
   - 4D max: no coalescence
   See also `../4_21_to_4d/README.md` and `../4_21_to_2d/README.md`.

