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
| Volume (unit roots) | 8.8800647944 |
| Volume (2× roots) | 142.0810367 |
| Hull vertices | 48 |
| Hull simplices (triangulated 3-faces) | 192 |
| Edges (nearest-neighbor) | 144 |
| Vertex degree | 6 (uniform) |
| Equilateral hull vertices | **No** (2 norm classes) |
| Generator coalescence | **None** (all 8 generators distinct) |
| **Symmetry group** | **W(F₄), order 1152** |

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

The 48 hull vertices split into **two orbits of 24** under W(F₄), each
forming a **24-cell** (the self-dual regular 4-polytope with 24 vertices,
96 edges, 96 triangular faces, 24 octahedral cells):

| | Orbit A (large) | Orbit B (small) |
|---|---|---|
| Vertices | 24 | 24 |
| |v|² | 1.911 | 1.411 |
| |v| | 1.383 | 1.188 |
| Within-orbit edges | 96 | 96 |
| Nearest neighbors | 8 per vertex | 8 per vertex |
| Normalized inner products | {−1, −½, 0, +½} | {−1, −½, 0, +½} |
| Matches 24-cell | **Yes** | **Yes** |
| E₈ integer roots on orbit | 8 | 12 |
| E₈ half-integer roots on orbit | 16 | 12 |

The two 24-cells are **dually positioned** (one is rotated relative to the
other), with the 144 nearest-neighbor edges of the full hull connecting
vertices **between** the two orbits (cross-orbit distance 1.0 < any
within-orbit distance). This is the same combinatorial structure as the
**F₄ root system** (48 roots = 24 long + 24 short, each orbit a 24-cell).

The norm² ratio is 1.911/1.411 ≈ 1.354, which differs from the standard
F₄ ratio of 2:1. This means the 4D max-volume projection produces a
non-standard-scale realization of the F₄ root arrangement, preserving the
combinatorial and group-theoretic structure but with a different metric.

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

## Files

- `../../max_shadow_421_to_4d.py` — main script. Builds the 240 E₈ roots,
  computes the H₄ baseline, runs basin-hopping over ℝ⁴ˣ⁸ matrices, and
  cross-checks against the old `e8.py` optimal bases.
- `../../e8.py` — the original Python 2 analysis (preserved for reference).
