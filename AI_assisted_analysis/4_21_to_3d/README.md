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

| | Max-volume (3D) | 7-fold (D₇d, 3D) | 4D max-volume (e8.py) | H₄ to 4D (600-cell) |
|---|---|---|---|---|
| Vertices | **18** (3 orbits of 6) | 240 | 48 | 600 |
| Edges    | **48** (4 orbits of 12) | (huge) | 144 | 720 |
| Facets   | **32** (all triangles) | (huge) | n/a | n/a |
| Symmetry | **D₃d, order 12** (with −I) | D₇d, order 28 | ? | H₄, order 14400 |
| Equilateral hull verts? | No (3 norm classes) | Yes | No | Yes |
| Volume (orthonormal frame) | 7.0553 | 6.2668 | 142.08 | 129.44 |

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
4. Does this 4-axis-coalescence pattern generalize? (E.g. for 4_21 → 4D,
   do some axes coalesce at the optimum? `e8.py` did not check.)
