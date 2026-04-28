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
| Equilateral hull vertices | **No** (2 norm classes) |
| Generator coalescence | **None** (all 8 generators distinct) |

Vertex norm classes:
- 24 vertices at |v|² ≈ 1.411 (|v| ≈ 1.188)
- 24 vertices at |v|² ≈ 1.911 (|v| ≈ 1.383)

Each hull vertex maps to a unique E₈ root (no projection degeneracy among
hull vertices).

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
| **Global max** | 8.8800648 | 48 | Non-equilateral |
| Second-best | 8.6942536 | — | About 2.1% below |
| H₄ baseline | 8.0901699 | 120 (600-cell) | 8.9% below |

The optimizer consistently finds the same max across 30 trials with
different random seeds. Trials that don't reach the global max converge
to the second-best at ≈ 8.694.

### Pairwise distance analysis (48 hull vertices)

| Distance | Count | Notes |
|----------|-------|-------|
| 1.0000 | 144 | Nearest-neighbor edges |
| 1.1880 | 96 | |
| 1.3825 | 96 | |
| 1.6801 | 72 | |
| 1.8229 | 288 | Most common |
| 1.9552 | 72 | |
| 2.0577 | 96 | |
| 2.3761 | 156 | |
| 2.3946 | 96 | |
| 2.7651 | 12 | Diameter (6 antipodal pairs) |

Total: C(48,2) = 1128 pairs. 10 distinct distance classes.

## Comparison across dimensions

| 4_21 → kD | Max volume | Hull vertices | Equilateral? | Coalescence |
|------------|-----------|---------------|--------------|-------------|
| 8 → 2 | 3√3 ≈ 5.196 | 6 | Yes | 1+5+1+1 (massive) |
| 8 → 3 | 8√7/3 ≈ 7.055 | 18 | No (3 norms) | 4 axes coalesce |
| **8 → 4** | **≈ 8.880** | **48** | **No (2 norms)** | **None** |

As dimension increases: hull vertex count grows (6 → 18 → 48), coalescence
decreases (massive → partial → none), and the hull shape becomes less
symmetric. The pattern that "optimal projections are less structured in
higher target dimensions" is consistent with the n-cube findings.

## Files

- `../../max_shadow_421_to_4d.py` — main script. Builds the 240 E₈ roots,
  computes the H₄ baseline, runs basin-hopping over ℝ⁴ˣ⁸ matrices, and
  cross-checks against the old `e8.py` optimal bases.
- `../../e8.py` — the original Python 2 analysis (preserved for reference).
