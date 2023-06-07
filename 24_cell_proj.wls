#!/usr/bin/env wolframscript

(* norm = Sqrt[3] 2 *)

(* a2 = { {1, 1, 2, 0}, {Sqrt[3], -Sqrt[3], 0, 0} };
f4 = { {1, 1, Sqrt[3] - 1, 0}, {Sqrt[3], -Sqrt[3], 0, Sqrt[3] - 3} };

Print[Norm[a2[[1]]]/Norm[f4[[1]] ] //FullSimplify ];
Print[Norm[a2[[1]]]/Norm[f4[[1]] ] //N ];

norm = Sqrt[2] Norm[f4[[1]] ]//FullSimplify;

a2 = a2/Norm[a2[[1]]]*norm//FullSimplify;
f4 = f4/Norm[f4[[1]]]*norm//FullSimplify;

Print[Map[Norm, a2 ]//N];
Print[Map[Norm, f4 ]//N];

Print[a2];
Print[f4];

Print[a2//N];
Print[f4//N];
 *)

proj = Transpose[{
	{Cos[t] a, Sin[t] a},
	{Cos[t + 120/180*Pi] a, Sin[t + 120/180*Pi] a},
	{Cos[t + 30/180*Pi] b, Sin[t + 30/180*Pi] b},
	{Cos[t + 150/180*Pi] b, Sin[t + 150/180*Pi] b}
}];

norms = Simplify[ Map[Norm[#]^2 &, proj], a > 0 && b > 0 && t > 0]/.{a -> 1, b -> 2};

Print[norms];

sol = Solve[norms[[1]] == norms[[2]], {t}];
Print[sol];


