# Maximum-Volume 4D Shadow of the 4_21 Polytope — Verification

## Setup

The **4_21 polytope** has the 240 E₈ roots as its vertices in ℝ⁸:
- **112 integer roots:** all signed permutations of (±1, ±1, 0, 0, 0, 0, 0, 0)
- **128 half-integer roots:** (±½, ±½, ±½, ±½, ±½, ±½, ±½, ±½) with an even
  number of minus signs

Every root has Euclidean norm √2.

We seek the 4-dimensional subspace W ⊂ ℝ⁸ that maximizes
**Vol₄(proj_W(4_21))** — i.e. the convex-hull volume of the 240 projected
points.

## Summary: verification of e8.py

The earlier `e8.py` (Python 2) found:
- **H₄ (phi-based) projection → 600-cell**: volume 129.44 (in the old 2×-scaled
  vertex convention)
- **Numerical maximum → 48-vertex polytope**: volume 142.08

Both are **confirmed** by our modern Python 3 re-verification, with a
correction to one error in the old comments.

### Scaling note

The old `e8.py` used vertices scaled by 2× (integer roots with ±2 instead
of ±1, half-integer roots with ±1 instead of ±½). All roots had norm 2√2
instead of √2. In 4D, this multiplies volumes by 2⁴ = 16:

| Quantity | Old e8.py (2× roots) | This work (unit roots) | Ratio |
|----------|---------------------|----------------------|-------|
| H₄ volume | 129.4427191 | 8.0901699437 | 16.00 |
| Max volume | 142.0810367 | 8.8800647944 | 16.00 |

## Detailed results

### H₄ (phi-based) projection: the 600-cell

| Property | Value |
|----------|-------|
| Volume (unit roots) | 8.0901699437 |
| Volume (2× roots) | 129.4427191 |
| Hull vertices | **120** (not 600) |
| Edges | 720 |
| Hull vertex norm | 1.203002 (uniform) |
| Equilateral | Yes |
| Symmetry | H₄ (order 14400) |

The projection basis from `e8.py` (φ = (1+√5)/2):

    basis₁ = (1, φ, 0, −1, φ, 0, 0, 0)
    basis₂ = (φ, 0, 1, φ, 0, −1, 0, 0)
    basis₃ = (0, 1, φ, 0, −1, φ, 0, 0)
    basis₄ = (0, 0, 0, 0, 0, 0, φ+1, φ−1)

### Max-volume 4D projection

| Property | Value |
|----------|-------|
| **Volume (unit roots)** | **(17 + 7√7)/4 ≈ 8.8800647944** |
| Volume (2× roots) | 68 + 28√7 ≈ 142.0810367 |
| Hull vertices | 48 |
| f-vector | (48, 240, 384, 192) |
| Edges | 240 = 144 (length 1, A–B) + 96 (length √((5+√7)/4), A–A) |
| Vertex degree | 20 (orbit A), 12 (orbit B) |
| Equilateral hull vertices | **No** (2 norm classes) |
| Generator coalescence | **None** (all 8 generators distinct) |
| **Symmetry group** | **W(F₄), order 1152** |

### Closed-form solution

The 48 hull vertices are the **F₄ root system** in a canonical coordinate frame:

**Orbit A (24 "long" roots):** all permutations of **(±a, ±a, 0, 0)** where

    a = √((5+√7)/8)

These form a **24-cell** of circumradius rA = √((5+√7)/4) ≈ 1.3825.

**Orbit B (24 "short" roots):** all permutations of **(±b, 0, 0, 0)** together
with all **(±b/2, ±b/2, ±b/2, ±b/2)** where

    b = √((3+√7)/4)

These form a **dual 24-cell** of circumradius rB = √((3+√7)/4) ≈ 1.1880.

**Key algebraic relations:**

    rA² = (5 + √7)/4
    rB² = (3 + √7)/4
    rA² − rB² = 1/2         (exact)
    rA²/rB² = 4 − √7        (norm² ratio)
    V = (17 + 7√7)/4         (volume, conjectured)

All values lie in Q(√7) — the same field as the 3D max volume 8√7/3.
The appearance of √7 in both the 3D and 4D optima suggests a deep connection
to the E₈ algebraic structure.

**Verification:** The standard F₄ root system constructed with these radii
reproduces all 10 pairwise distance classes of the numerical hull exactly,
including the cross-orbit minimum distance of 1.0. The volume matches the
conjectured (17+7√7)/4 to 13 significant digits.

**Orbit comparison:**

| | Orbit A ("long") | Orbit B ("short") |
|---|---|---|
| **Vertex form** | perms of (±a, ±a, 0, 0) | perms of (±b, 0, 0, 0) + (±b/2)⁴ |
| \|v\|² | (5+√7)/4 ≈ 1.911 | (3+√7)/4 ≈ 1.411 |
| # vertices | 24 | 24 |
| Within-orbit edges | 96 | 96 |
| Within-orbit distances | 4 classes | 4 classes |
| E₈ integer roots | 8 | 12 |
| E₈ half-integer roots | 16 | 12 |

**Orientation.** The two 24-cells are in the standard F₄ **dual**
orientation: orbit A consists of "D₄ roots" (perms of (±a,±a,0,0)) and
orbit B consists of "dual D₄ roots" (axis vertices (±b,0,0,0) plus all-half
vertices (±b/2)⁴). The relative rotation between them is a 45° rotation in
a 2D plane combined with a reflection — the D₄ triality transformation.

**Non-standard scale.** The norm² ratio rA²/rB² = 4 − √7 ≈ 1.354 differs
from both the standard F₄ ratio (2:1) and the equal-radii ratio (1:1). The
4D max-volume projection produces a non-standard-scale realization of the
F₄ root system, with the two 24-cell radii determined by the optimization
and lying in Q(√7). The 144 nearest-neighbor (length-1) edges connect
vertices **between** the two orbits — the same combinatorial pattern as the
F₄ root system (48 roots = 24 long + 24 short, each orbit a 24-cell).

### Explicit 4×8 projection basis

A canonical 4×8 orthonormal projection M whose action on the 240 E₈ roots
yields the F₄ hull above (in the canonical 4D frame) has all entries equal
to ±c or ±d with

    c = √((5 + √7)/32) = a/2          ≈ 0.488805
    d = √((3 − √7)/32) = 1/(8b)        ≈ 0.105215

(equivalently c² + d² = 1/4 and c² − d² = (1+√7)/16). One such matrix:

    basis₁ = ( +d, −c, −d, +c, −d, −d, −c, −c )
    basis₂ = ( +d, −c, +d, +c, −d, +d, +c, +c )
    basis₃ = ( −c, +d, +c, +d, −c, −c, +d, −d )
    basis₄ = ( +d, −c, +d, −c, +d, −d, +c, −c )

Properties:
- Each basis vector has 4 entries of magnitude c and 4 of magnitude d;
  norm² = 4c² + 4d² = 1 (the four basis vectors are orthonormal, M Mᵀ = I).
- Each column has either {1·c, 3·d} (columns 1,3,5,6) — call these "small"
  generators with norm² = c² + 3d² = (7−√7)/16 ≈ 0.272, or
  {3·c, 1·d} (columns 2,4,7,8) — "large" generators with norm² =
  3c² + d² = (9+√7)/16 ≈ 0.728.
- Total Σ_j |g_j|² = 4·(7−√7)/16 + 4·(9+√7)/16 = 64/16 = 4 ✓
  (= ambient dim 4 since the basis is orthonormal).
- The 4+4 split into two generator-norm classes matches the "no generator
  coalescence but two norm classes" observation noted earlier.

The matrix is unique only up to W(F₄) (order 1152) acting on the basis
vectors (rotations/reflections in the 4D image) and S₈ × {±1}⁸ on columns
(relabelling/sign-flipping of E₈ axes). Many sign patterns give equally
valid bases — the table above is one specific representative.

### Symmetry: the Weyl group of F₄

The symmetry group of the 48-vertex hull is **W(F₄)**, the Weyl group of the
exceptional Lie algebra F₄, with order **1152 = 2⁷ · 3²**. This was found by
brute-force enumeration of all orthogonal matrices in O(4) preserving the
vertex set. The group decomposes as

    W(F₄) = G⁺ × {I, −I}

where G⁺ (the rotation subgroup) has order 576.

Element orders and counts:

| Order | Count |
|-------|-------|
| 1 | 1 |
| 2 | 139 |
| 3 | 80 |
| 4 | 228 |
| 6 | 464 |
| 8 | 144 |
| 12 | 96 |
| **Total** | **1152** |

### f-vector and combinatorics

For the max-volume hull (rA²/rB² = 4 − √7):

| k | element | count |
|---|---|---|
| 0 | vertices         | **48**  (24 in orbit A + 24 in orbit B) |
| 1 | edges            | **240** |
| 2 | 2-faces (triangles) | **384** |
| 3 | 3-cells (tetrahedra)| **192** |

Euler check: V − E + F − C = 48 − 240 + 384 − 192 = 0 ✓.

**Edges (2 length classes):**

| length (exact) | length (≈) | count | type |
|---|---|---|---|
| **1** | 1.0000 | 144 | A–B (cross-orbit, nearest-neighbor) |
| **√((5+√7)/4)** | 1.3825 | 96 | A–A (24-cell edges of orbit A) |

No B–B edges. Every orbit-A vertex has 12 A–B edges + 8 A–A edges = degree 20.
Every orbit-B vertex has 12 A–B edges + 0 within = degree 12.

**2-faces (all triangles, 2 shapes):**

| shape | side lengths | count |
|---|---|---|
| isoceles | (1, 1, √((5+√7)/4))  | 288 |
| equilateral | (√((5+√7)/4))³    | 96  |

The 96 equilateral triangles are exactly the 96 triangular faces of the
orbit-A 24-cell.

**3-cells (all 192 are tetrahedra of type 3A·1B):**

Every cell has 3 orbit-A vertices forming an equilateral triangle (long
edges) plus 1 orbit-B apex connected to the base by 3 unit-length A–B edges.
This is an isoceles ("pyramidal") tetrahedron, not a regular one. Each
orbit-B vertex is the apex of 192/24 = 8 cells; each orbit-A vertex is in
192·3/24 = 24 cells.

Geometrically: take the orbit-A 24-cell and erect a tetrahedral cap on
**both sides** of each of its 96 triangular faces (96 × 2 = 192 caps). Each
cap's apex is an orbit-B vertex. The orbit-B 24-cell sits inside but its
own faces and edges are not visible on the boundary.

### Cell shape: right regular trigonal pyramid

All 192 cells are **congruent** (same edge multiset: three of length 1 and
three of length L = √((5+√7)/4) ≈ 1.3825). Each is a **right pyramid over
an equilateral triangle** — i.e. the orbit-B apex sits exactly above the
centroid of the equilateral A-triangle base. Symmetry group of one cell is
**C₃ᵥ** (order 6).

| Quantity | Closed form | Value |
|---|---|---|
| Base edge (A–A) | L = √((5+√7)/4) | 1.3825 |
| Apex edge (A–B) | 1 (exactly) | 1.0000 |
| Apex height above base | h = √((7−√7)/12) | 0.6024 |
| Apex angle (at B, between two apex edges) | arccos((3−√7)/8) | 87.462° |
| Vertex angle at A in side face | arccos(L/2) | 46.269° |
| Dihedral along A–A edge | arccos(0.5523) | 56.473° |
| Dihedral along A–B edge | arccos(0.0424) | 87.570° |
| Cell 3-volume | V_cell = √((77+19√7)/4608) | 0.1662 |

The apex angle and the A–B dihedral are both close to (but not exactly) 90°
— both would equal 90° if L² = 2; here L² = (5+√7)/4 ≈ 1.911. The cell is
**not a regular tetrahedron**, **not a classical isoceles tetrahedron**
(opposite-edge-equal), and **not orthocentric** — just a generic regular
trigonal pyramid scaled by the volume-maximizing edge ratio. Hull volume
check: with all 192 cells congruent and W(F₄)-equivalent, the inradius
(common distance from origin to each cell hyperplane) is r_in² = (35+13√7)/56,
and (1/4)·192·V_cell·r_in = (17+7√7)/4 ✓.

### Comparison with the disphenoidal 288-cell

The **disphenoidal 288-cell** is the convex hull of the same 48-vertex
union of two dual 24-cells when **both 24-cells share a common circumradius**
(rA = rB). It is a uniform polytope, vertex-transitive under W(F₄), and is
the dual of the bitruncated 24-cell.

| Feature | Disphenoidal 288-cell (rA = rB) | Max-volume hull (rA² = (4−√7)·rB²) |
|---|---|---|
| f-vector | (48, **336**, 576, 288)        | (48, **240**, 384, 192) |
| edges    | 144 short + **192 long** = 336 | 144 short + **96 long** = 240 |
| edge ratio long/short | 1.000 / 0.7654 ≈ 1.307 | 1.383 / 1.000 = 1.383 |
| 2-face shapes | 576 isoceles (one shape)   | 288 isoceles + 96 equilateral |
| 3-cells  | 288 disphenoids (irregular tets) | 192 isoceles tetrahedra |
| Symmetry on vertices | transitive (one orbit of 48) | two orbits of 24 |
| Cells per A-vertex / B-vertex | 24 / 24 (transitive) | 24 / 8 |
| Field of coordinates | Q (rational up to scale) | Q(√7) |

Why the long-edge count drops from 192 → 96: when rA = rB, every "long"
edge of the disphenoidal 288-cell connects two A-vertices that are 24-cell
edge-neighbors **across** an A–A "diagonal pair" of the outer 24-cell. As
rA increases past rB, the orbit-A vertices push outward; pairs that were
just barely connected on the rA = rB hull now have midpoints lying interior
to the orbit-B vertex on the same radial direction, and they cease to be
edges. Specifically, the midpoint of the would-be edge from (a,a,0,0) to
(a,−a,0,0) is (a,0,0,0); it lies inside the hull iff a < b, i.e. iff
rA²/2 < rB², i.e. iff rA²/rB² < 2. Since 4 − √7 ≈ 1.354 < 2, all 96
"face-diagonal" A–A pairs in the outer 24-cell drop off the hull.

The disphenoidal 288-cell has higher symmetry and a tidier combinatorial
structure (fully transitive on vertices, edges, and cells), but **smaller
projected volume**: it corresponds to an F₄-symmetric subspace where the
two 24-cell scales are forced equal, which is not the volume-maximizing
configuration. The optimization breaks the rA ↔ rB symmetry and pays for
the loss of uniformity with extra volume.

### Why F₄ appears: E₈ ⊃ F₄

The appearance of F₄ is deeply connected to the structure of E₈.
The E₈ root system contains F₄ as a maximal root subsystem: there exist
48 E₈ roots forming an F₄ root system in a 4-dimensional subspace of ℝ⁸.
The max-volume 4D projection appears to be related to (though not identical
with) such an F₄ subsystem direction. The 48 hull vertices are 48 specific
E₈ roots that project to the boundary, arranged as two 24-cells.

### Generator structure

The 4+4 split of the 8 generators into two norm classes (|g| = √((9+√7)/16) ≈
0.853 and √((7−√7)/16) ≈ 0.522, see the explicit basis above) is reminiscent
of the block decomposition V = V₁ ⊕ V₂ found in the 8→4 max-shadow
**zonotope** (8-cube projection), though the 4_21 is not a zonotope so the
structural interpretation differs.

### Local maxima

| Local max | Volume (unit roots) | Hull vertices | Notes |
|-----------|-------------------|---------------|-------|
| **Global max** | 8.8800648 | 48 | F₄ arrangement, W(F₄) symmetry |
| Second-best | 8.6942536 | — | About 2.1% below |
| H₄ baseline | 8.0901699 | 120 (600-cell) | 8.9% below |

The optimizer consistently finds the same max across 30 trials with
different random seeds. Trials that don't reach the global max converge
to the second-best at ≈ 8.694.

### Pairwise distance analysis (48 hull vertices)

| Distance | Count | Notes |
|----------|-------|-------|
| 1.0000 | 144 | Edges (all cross-orbit) |
| 1.1880 | 96 | Within orbit B |
| 1.3825 | 96 | Within orbit A |
| 1.6801 | 72 | Within orbit B |
| 1.8229 | 288 | Cross-orbit |
| 1.9552 | 72 | Within orbit A |
| 2.0577 | 96 | Within orbit B |
| 2.3761 | 156 | Cross-orbit + within B |
| 2.3946 | 96 | Within orbit A |
| 2.7651 | 12 | Diameter (within orbit A) |

Total: C(48,2) = 1128 pairs. 10 distinct distance classes.

## Subspace nesting

The optimal subspaces for different target dimensions do **not** nest
monotonically:

    2D (A₂ plane) ⊂ 3D (D₃d subspace) ⊄ 4D (F₄ subspace)

- The 2D max (3√3) is achieved within the 3D max subspace (the xy-plane
  of the canonicalized 3D projection is an A₂ plane).
- The 3D max (8√7/3) is **not** achievable within the 4D max subspace:
  the best 3D volume within the 4D space is 7.033 (99.7% of the true max).
- The 2D max is also not achievable within the 4D subspace: best 2D area
  within 4D is 4.966 (95.6% of 3√3).

## Comparison across dimensions

| 4_21 → kD | Max volume | Hull vertices | Symmetry | Coalescence |
|------------|-----------|---------------|----------|-------------|
| 8 → 2 | 3√3 ≈ 5.196 | 6 | D₆, order 12 | 1+5+1+1 (massive) |
| 8 → 3 | 8√7/3 ≈ 7.055 | 18 | D₃d, order 12 | 4 axes coalesce |
| **8 → 4** | **≈ 8.880** | **48** | **W(F₄), order 1152** | **None** |

Symmetry jumps dramatically from 3D (order 12) to 4D (order 1152). The 4D
hull is far more symmetric than the 3D hull despite having more vertices.
The F₄ connection suggests the 4D optimum is governed by the E₈ ⊃ F₄
algebraic structure rather than the lower-dimensional geometric constraints
that shape the 2D and 3D optima.

## vZome impossibility

Like the 3D max-volume hull, the 4D F₄ hull is **not constructible in any of
vZome's standard fields**. The 48 vertex coordinates depend on

- a² = (5 + √7)/8   (orbit A)
- b² = (3 + √7)/4   (orbit B)

and the ratio (b/a)² = (8 + 2√7)/9 still lies in Q(√7), so no global rescaling
can clear √7. Projecting to 3D via either the long-root or short-root view
introduces only √2, √3, √6 from the basis, leaving the projected coordinates
in Q(√2, √3, √7).

vZome's built-in fields — rational, golden Q(φ), root2 Q(√2), root3 Q(√3),
heptagonal Q(2cos π/7), snubDodec, snubCube — none contain √7. (The
heptagonal field is a totally real cubic with discriminant 49, but its only
quadratic subfield is Q(√−7), not Q(√7).) Since the field requirement is
isometry-invariant, no rotation or choice of projection direction can remove
√7. The smallest field that supports both the 3D and 4D max-volume hulls is
Q(√2, √3, √7), which would have to be added to vZome as a custom field.

## Files

- `../../max_shadow_421_to_4d.py` — main script. Builds the 240 E₈ roots,
  computes the H₄ baseline, runs basin-hopping over ℝ⁴ˣ⁸ matrices, and
  cross-checks against the old `e8.py` optimal bases.
- `../../e8.py` — the original Python 2 analysis (preserved for reference).
- `generate_stl.py` — builds ball-and-stick STL wireframes from the closed-form
  F₄ vertices, projected to 3D via two dual viewpoints.
- `max_volume_421_to_4d_short_root_view.stl` — 4D→3D projection along the
  short-root direction (1,1,1,1)/2. An orbit B vertex at the pole. S₄ symmetry
  (order 24) in the 3D projection.
- `max_volume_421_to_4d_long_root_view.stl` — 4D→3D projection along the
  long-root direction (1,1,0,0)/√2. An orbit A vertex at the pole. D₄×Z₂
  symmetry (order 16) in the 3D projection. This is the "dual" viewpoint
  where the roles of the two 24-cells are exchanged.

### Projection coordinates (4D → 3D, for STL views)

The two STL views use the following orthonormal 3D bases within the 4D
space (columns of the 4×3 projection matrix):

**Short-root view** (project along (1,1,1,1)/2):

    x = (1, −1, 0, 0) / √2
    y = (1, 1, −2, 0) / √6
    z = (1, 1, 1, −3) / √12

**Long-root view** (project along (1,1,0,0)/√2):

    x = (1, −1, 0, 0) / √2
    y = (0, 0, 1, 0)
    z = (0, 0, 0, 1)

Both views produce 7 depth layers with 33 distinct 3D points. In the
short-root view, orbit B vertices sit at the deepest layers (pole/antipole);
in the long-root view, orbit A vertices take that role — a clean swap
reflecting the F₄ long/short root duality.
