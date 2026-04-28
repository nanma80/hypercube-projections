# Maximum-Volume 2D Shadow of the 4_21 Polytope

## Setup

The **4_21 polytope** has the 240 E₈ roots as its vertices in ℝ⁸:
- **112 integer roots:** all signed permutations of (±1, ±1, 0, 0, 0, 0, 0, 0)
- **128 half-integer roots:** (±½, ±½, ±½, ±½, ±½, ±½, ±½, ±½) with an even
  number of minus signs

Every root has Euclidean norm √2.

We seek the 2-dimensional subspace W ⊂ ℝ⁸ that maximizes
**Area(proj_W(4_21))** — i.e. the convex-hull area of the 240 projected
points.

## Result

**Max area = 3√3 ≈ 5.1961524227**, achieved by projecting onto any
**A₂ root subsystem plane** of E₈. The convex hull is a **regular hexagon**
inscribed in a circle of radius √2.

Verified by 40 basin-hopping trials, all converging to either 3√3 (the
global max) or the local max 2 + 2√2 ≈ 4.8284 (the 8-cube's max 2D
shadow).

| Projection                                | Area      | Ratio  |
|-------------------------------------------|-----------|--------|
| **Max-area (A₂ plane, this work)**        | 5.1962    | 1.000  |
| B₈ Petrie plane (regular 8-gon)           | 4.6447    | 0.894  |
| 2 + 2√2 (8-cube max, local max for 4_21) | 4.8284    | 0.929  |
| D₈ Petrie plane                           | 4.5306    | 0.872  |
| E₈ phi-based 2D (first 2 rows of H₄)     | 4.0652    | 0.782  |

### Hull combinatorics

| Property | Value |
|----------|-------|
| Area | 3√3 = 5.1961524227 |
| f-vector | (6, 6) |
| Shape | **Regular hexagon** |
| Circumradius | √2 |
| Edge length | √2 |
| Symmetry | D₆ (dihedral, order 12) |
| Hull vertices equilateral | Yes (all at norm √2) |

## The optimal 2D plane is an A₂ root subsystem plane

The E₈ root system contains **A₂** as a root subsystem: any three E₈ roots
{r₁, r₂, r₃} at mutual 60° angles (inner product 1) with r₃ = r₁ − r₂
generate a 2D plane containing 6 roots forming a regular hexagon:
{±r₁, ±r₂, ±r₃}.

For example, r₁ = e₁ + e₂ = (1,1,0,0,0,0,0,0), r₂ = e₁ − e₃ = (1,0,−1,0,0,0,0,0),
r₃ = e₂ + e₃ = (0,1,1,0,0,0,0,0). These are all E₈ integer roots at mutual
inner product 1 (= |r|²cos 60° = 2 · ½).

Projecting all 240 E₈ roots onto this 2D plane:
- The 6 A₂ roots project to the hexagon vertices (at full norm √2)
- All other 234 roots project **strictly inside** the hexagon
- 240 roots collapse to only **13 distinct 2D points**
- **72 roots project to the origin** (those orthogonal to the A₂ plane)

This is the maximum because:
1. The hexagon circumradius √2 equals the E₈ root norm — no projection can
   have vertices farther from the origin
2. A regular hexagon maximizes area among all convex sets inscribed in a
   circle (for a fixed number of vertices ≤ 6), and no 2D projection of
   the E₈ roots places more than 6 points at the maximum radius

All A₂ subsystem planes of E₈ are conjugate under the Weyl group W(E₈),
so the maximum is achieved by any such plane and the orbit under W(E₈) gives
all optimal projections.

## Generator coalescence: 1 + 5 + 1 + 1 pattern

The optimal 2×8 projection matrix Q has generators g_i = Q[i,:] (row i)
falling into just 2 directions:

    g₀ = +a                          |a| = 1/√6
    g₁ = g₄ = g₅ = g₆ = g₇ = −a    (5 coalesced axes)
    g₂ = +b                          |b| = 1/√2
    g₃ = −b

where a ⊥ b. The kernel of the projection is 6-dimensional (rank 2 from
rank 8), with a highly symmetric 4-dimensional subspace arising from the
5-fold coalescence of axes 1,4,5,6,7.

## The 2D max-area plane nests inside the 3D max-volume subspace

The 3D max-volume shadow (see `../4_21_to_3d/README.md`) has an equatorial
regular hexagon at z = 0 with 6 vertices at norm √2. Projecting the 3D
shadow to its xy-plane gives **exactly the same regular hexagon** with area
3√3 — the 2D maximum.

This works because the xy-plane of the canonicalized 3D projection IS an A₂
root subsystem plane. The 8 generators projected to xy are:

    g₀_xy = (+1/√2, 0),  g₁_xy = (−1/√2, 0)       [|g| = 1/√2]
    g₂_xy = (0, −1/√6)                              [|g| = 1/√6]
    g₃_xy = g₄_xy = g₅_xy = g₆_xy = g₇_xy = (0, +1/√6)  [|g| = 1/√6]

The two E₈ roots ±(e₁ − e₂) project to (±√2, 0) on the hexagon, and four
half-integer roots project to the remaining four hexagonal vertices at 60°
intervals — confirming A₂ structure.

**The nesting is one-way.** The 2D max embeds in the 3D max, but the 3D
max does NOT embed in the 4D max: the best 3D volume achievable within
the 4D max-volume subspace is only 7.033 (99.7% of the true 3D max 7.055).
Similarly, the best 2D area within the 4D subspace is 4.966 (95.6% of 3√3).
The three optimal subspaces are:

    2D (A₂ plane) ⊂ 3D (D₃d subspace) ⊄ 4D (F₄ subspace)

## Comparison across dimensions

| 4_21 → kD | Volume/Area | Hull vertices | Shape | Symmetry |
|------------|-------------|---------------|-------|----------|
| **8 → 2** | **3√3 ≈ 5.196** | **6** | **Regular hexagon** | **D₆, order 12** |
| 8 → 3 | 8√7/3 ≈ 7.055 | 18 | Tented squashed cuboctahedron | D₃d, order 12 |
| 8 → 4 | (17+7√7)/4 ≈ 8.880 | 48 | Two 24-cells (F₄ root arrangement) | W(F₄), order 1152 |

## Files

- `../../max_shadow_421_to_2d.py` — main script. Builds the 240 E₈ roots,
  computes baselines (B₈ Petrie, D₈ Petrie, phi-based), and runs
  basin-hopping over ℝ²ˣ⁸ matrices to find the max-area 2D shadow.
