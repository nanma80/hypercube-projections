(* 6->3 Rhombic Triacontahedron: symbolic volume *)
(* 6 generators = 5-fold axes of icosahedron, unit vectors *)
(* These form an equilateral tight frame with PP^T = 2I *)

phi = (1 + Sqrt[5])/2;
N0 = Sqrt[(21 + 9 Sqrt[5])/2];

(* Generators: half the icosahedron vertices, normalized *)
(* Use the standard icosahedron vertex set *)
gens = {
  {0, 1, phi}, {0, 1, -phi}, {1, phi, 0},
  {1, -phi, 0}, {phi, 0, 1}, {phi, 0, -1}
} / Sqrt[1 + phi^2];

Print["Generator norms: ", FullSimplify[Map[Norm, gens]]];
Print["PP^T = ", FullSimplify[Sum[Outer[Times, gens[[i]], gens[[i]]], {i, 1, 6}]]];

(* Shadow volume = sum |det(P_S)| over C(6,3)=20 triples *)
vol = Sum[
  Abs[Det[{gens[[i]], gens[[j]], gens[[k]]}]],
  {i, 1, 6}, {j, i+1, 6}, {k, j+1, 6}
];

volSimp = FullSimplify[vol, Assumptions -> {Sqrt[5] > 0}];
Print["V = ", volSimp];
Print["V (numeric) = ", N[volSimp, 20]];
Print["Expected: 4.3525017990"];
