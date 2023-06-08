#!/usr/bin/env wolframscript

epsilon = 0.0000000001;

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

(* proj[t_, a_, b_, ra_, rb_] := Transpose[{
				{Cos[t] a, Sin[t] a},
				{Cos[t + k1/180*Pi] a ra, Sin[t + k1/180*Pi] a ra},
				{Cos[t + k3/180*Pi] b, Sin[t + k3/180*Pi] b},
				{Cos[t + (k3+k2)/180*Pi] b rb, Sin[t + (k3+k2)/180*Pi] b rb}
			}]//N;
 *)

For[k1 = 120, k1 <= 360, k1 += 60,
	If[Mod[k1, 180] == 0, Continue[]];
	For[k2 = 60, k2 <= 360, k2 += 60,
		If[Mod[k2, 180] == 0, Continue[]];
		For[k3 = 30, k3 <= 330, k3 += 60,
			(* If[{k1, k2, k3} != {240, 120, 150}, Continue[]]; *)
			proj[t_, a_, b_, ra_, rb1_, rb2_] := Transpose[{
				{Cos[t] a, Sin[t] a},
				{Cos[t + k1/180*Pi] a ra, Sin[t + k1/180*Pi] a ra},
				{Cos[t + k3/180*Pi] b rb1, Sin[t + k3/180*Pi] b rb1},
				{Cos[t + (k3+k2)/180*Pi] b rb2, Sin[t + (k3+k2)/180*Pi] b rb2}
			}]//N;

			aValue = 1;
			bValue = 1.1;
			(* bValue = Sqrt[3]; *)
			(* ratio = 1/Sqrt[2]; *)

(* 			sol = FindRoot[
				{Dot[proj[t, aValue, bValue, ra, rb1, rb2][[1]], proj[t, aValue, bValue, ra, rb1, rb2][[1]] ] == Dot[proj[t, aValue, bValue, ra, rb1, rb2][[2]], proj[t, aValue, bValue, ra, rb1, rb2][[2]] ],
				Dot[proj[t, aValue, bValue, ra, rb1, rb2][[1]], proj[t, aValue, bValue, ra, rb1, rb2][[2]] ] == 0}
			, {{t, 0}, {ratio, 1}}
			, AccuracyGoal -> 4, PrecisionGoal -> 4]; *)

			sol = FindInstance[
				Dot[proj[t, aValue, bValue, ra, rb1, rb2][[1]], proj[t, aValue, bValue, ra, rb1, rb2][[1]] ] == Dot[proj[t, aValue, bValue, ra, rb1, rb2][[2]], proj[t, aValue, bValue, ra, rb1, rb2][[2]] ] &&
				Dot[proj[t, aValue, bValue, ra, rb1, rb2][[1]], proj[t, aValue, bValue, ra, rb1, rb2][[2]] ] == 0
				&& ra > 0 && rb1 > 0 && rb2 > 0 && t > 0 && t < 2 Pi
			, {t, ra, rb1, rb2}, Reals];

			If[Length[sol] != 1, Continue[]];

			Print[];

			Print[{{k1, k2, k3}, {0, k1, k3, k2+k3}}];

			(* Print[sol]; *)
			Print[sol//N];
			Print[N[{Dot[proj[t, aValue, bValue, ra, rb1, rb2][[1]], proj[t, aValue, bValue, ra, rb1, rb2][[1]] ], Dot[proj[t, aValue, bValue, ra, rb1, rb2][[2]], proj[t, aValue, bValue, ra, rb1, rb2][[2]] ],
				Dot[proj[t, aValue, bValue, ra, rb1, rb2][[1]], proj[t, aValue, bValue, ra, rb1, rb2][[2]] ]}/.sol]];
			(* If[Length[sol] == 0, Continue[]]; *)

			Print[N[{t/Pi, ra, rb1, rb2}/.sol]];
		];
	];
];

(* proj[t_, a_, b_] := Transpose[{
	{Cos[t] a, Sin[t] a},
	{Cos[t + 120/180*Pi] a, Sin[t + 120/180*Pi] a},
	{Cos[t + 30/180*Pi] b, Sin[t + 30/180*Pi] b},
	{Cos[t + (30+60 * 4)/180*Pi] b, Sin[t + (30+60 * 4)/180*Pi] b}
}]//N;

(* projection[a_, b_] := proj[-Pi/12, a, b]*Sqrt[8]; *)

aValue = 1;
bValue = 2;

sol = Solve[
	Dot[proj[t, aValue, bValue][[1]], proj[t, aValue, bValue][[1]] ] == Dot[proj[t, aValue, bValue][[2]], proj[t, aValue, bValue][[2]] ] &&
	Dot[proj[t, aValue, bValue][[1]], proj[t, aValue, bValue][[2]] ] == 0
, {t}];

Print[N[t/Pi]/.sol];
 *)



