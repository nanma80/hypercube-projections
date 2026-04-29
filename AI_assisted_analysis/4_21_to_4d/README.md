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

### Correction: old e8.py comment error

The old `e8.py` comments stated:

> "convex hull of the known bases: 600 cell. **600 vertices**, 720 edges"

This is **incorrect**. The 600-cell has **120 vertices**, 720 edges, 1200
triangular faces, and 600 tetrahedral cells. The "600" in "600-cell" refers
to the number of cells, not vertices. Our verification confirms:

- **120 hull vertices** (not 600)
- 720 edges ✓
- All hull vertices equilateral (norm 1.203 in orthonormal frame) ✓
- Only 120 of the 240 E₈ roots lie on the hull boundary; the other 120
  project to interior points of the 600-cell

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
| Hull simplices (triangulated 3-faces) | 192 |
| Edges (nearest-neighbor) | 144 |
| Vertex degree | 6 (uniform) |
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

### The two 24-cells: sizes and orientation

The two 24-cells are **different sizes** (rA/rB ≈ 1.164, as opposed to √2 ≈ 1.414
in the standard F₄ root system). They are in the **standard F₄ dual orientation**:
orbit A consists of "D₄ roots" (permutations of (±a,±a,0,0)) and orbit B consists
of "dual D₄ roots" (axis vertices (±b,0,0,0) plus all-half vertices (±b/2)⁴).
The relative rotation between the two 24-cells involves a 45° rotation in a
2D plane combined with a reflection, which is exactly the D₄ triality
transformation relating the D₄ root system to its dual.

| | Orbit A ("long") | Orbit B ("short") |
|---|---|---|
| **Vertex form** | perms of (±a, ±a, 0, 0) | perms of (±b, 0, 0, 0) + (±b/2)⁴ |
| |v|² | (5+√7)/4 ≈ 1.911 | (3+√7)/4 ≈ 1.411 |
| # vertices | 24 | 24 |
| Within-orbit edges | 96 | 96 |
| Within-orbit distances | 4 classes | 4 classes |
| E₈ integer roots | 8 | 12 |
| E₈ half-integer roots | 16 | 12 |

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

### The hull is two 24-cells: an F₄ root arrangement

The 48 hull vertices split into **two orbits of 24** under W(F₄). The 144
nearest-neighbor edges of the full hull connect vertices **between** the two
orbits (cross-orbit distance 1.0 < any within-orbit distance). This is the
same combinatorial structure as the **F₄ root system** (48 roots = 24 long +
24 short, each orbit a 24-cell).

The norm² ratio rA²/rB² = 4−√7 ≈ 1.354 differs from the standard F₄ ratio
of 2:1. The 4D max-volume projection produces a non-standard-scale realization
of the F₄ root system, where the two 24-cell radii are determined by the
optimization and lie in Q(√7).

### Why F₄ appears: E₈ ⊃ F₄

The appearance of F₄ is deeply connected to the structure of E₈.
The E₈ root system contains F₄ as a maximal root subsystem: there exist
48 E₈ roots forming an F₄ root system in a 4-dimensional subspace of ℝ⁸.
The max-volume 4D projection appears to be related to (though not identical
with) such an F₄ subsystem direction. The 48 hull vertices are 48 specific
E₈ roots that project to the boundary, arranged as two 24-cells.

### Generator structure

The 8 generators fall into **2 norm classes** (but no coalescence):
- 4 generators with |g| ≈ 0.853 (axes 0, 4, 6, 7 in one run)
- 4 generators with |g| ≈ 0.522 (axes 1, 2, 3, 5)

This 4+4 split is reminiscent of the block decomposition V = V₁ ⊕ V₂ found
in the 8→4 max-shadow **zonotope** (8-cube projection), though the 4_21 is
not a zonotope so the structural interpretation differs.

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

### Projection coordinates

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
