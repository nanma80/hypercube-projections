(* ================================================================
   SYMBOLIC VERIFICATION: 8->3 max-volume zonotope
   
   Derive V(h) exactly, differentiate, solve dV/dh=0,
   and verify h^2 = (2/21)(1 - rho + rho^2) where rho = 2cos(pi/7).
   ================================================================ *)

Print["================================================================"];
Print["SYMBOLIC DERIVATION OF V(h) AND OPTIMAL h"];
Print["================================================================"];

(* ---- Step 1: Volume formula ---- *)
(* 
  Generators:
    g_k = (r cos(2pi k/7), r sin(2pi k/7), h)  for k=0,...,6
    g_7 = (0, 0, b)
  where r = 4/Sqrt[21], b = Sqrt[8/3 - 7h^2].

  Shadow volume = (1/(8/3)^(3/2)) * Sum |det(P_S)| over 3-subsets S.
  
  Case 1 (all cone, i<j<k<7):
    det([g_i | g_j | g_k]) = r^2 h (sin(d_ij) + sin(d_jk) - sin(d_ik))
    where d_ij = 2pi(j-i)/7, d_jk = 2pi(k-j)/7, d_ik = 2pi(k-i)/7.
    
  Case 2 (two cone + axial, i<j<7):
    det([g_i | g_j | g_7]) = r^2 b sin(d_ij)
    where d_ij = 2pi(j-i)/7.
*)

(* Compute S1 = sum of |sin(d_ij)+sin(d_jk)-sin(d_ik)| over 35 triples *)
S1exact = Sum[
  Abs[Sin[2 Pi (j-i)/7] + Sin[2 Pi (k-j)/7] - Sin[2 Pi (k-i)/7]],
  {i, 0, 6}, {j, i+1, 6}, {k, j+1, 6}
];

(* Compute S2 = sum of |sin(2pi(j-i)/7)| over 21 pairs *)
S2exact = Sum[
  Abs[Sin[2 Pi (j-i)/7]],
  {i, 0, 6}, {j, i+1, 6}
];

Print["\nStep 1: Trigonometric sums"];
Print["  S1 (35 triple terms) = ", FullSimplify[S1exact]];
Print["  S2 (21 pair terms)   = ", FullSimplify[S2exact]];
Print["  S1 (numeric) = ", N[S1exact, 30]];
Print["  S2 (numeric) = ", N[S2exact, 30]];

(* ---- Step 2: Volume as function of h ---- *)
(* V_zonotope(h) = r^2 h S1 + r^2 b S2  where b = Sqrt[8/3 - 7h^2] *)
(* Shadow volume = V_zonotope / (8/3)^(3/2) *)
(* Since the normalization factor is constant, optimizing V_zonotope suffices. *)

Print["\nStep 2: Volume function"];
Print["  V_zono(h) = r^2 * (h * S1 + Sqrt[8/3 - 7h^2] * S2)"];
Print["  Shadow volume = V_zono / (8/3)^(3/2)"];

(* ---- Step 3: Differentiate and solve ---- *)
(* dV/dh = r^2 * (S1 + S2 * d/dh[Sqrt[8/3 - 7h^2]])
         = r^2 * (S1 - 7h*S2/Sqrt[8/3 - 7h^2])
         = r^2 * (S1 - 7h*S2/b) *)
(* Setting dV/dh = 0:
     S1 = 7h*S2/b
     S1*b = 7h*S2
     S1^2*b^2 = 49h^2*S2^2
     S1^2*(8/3 - 7h^2) = 49h^2*S2^2
     8*S1^2/3 = h^2*(7*S1^2 + 49*S2^2)
     h^2 = 8*S1^2 / (3*(7*S1^2 + 49*S2^2))
*)

hSquared = 8 S1exact^2 / (3 (7 S1exact^2 + 49 S2exact^2));

Print["\nStep 3: Solve dV/dh = 0"];
Print["  h^2 = 8*S1^2 / (3*(7*S1^2 + 49*S2^2))"];
Print["  h^2 (numeric) = ", N[hSquared, 30]];
Print["  h   (numeric) = ", N[Sqrt[hSquared], 30]];

(* ---- Step 4: Simplify using rho = 2cos(pi/7) ---- *)
(* S1 and S2 involve sin(k*pi/7) terms. These are algebraic over Q(rho). *)
(* sin(pi/7) = Sqrt[(2+rho)/4] * ... but let's let Mathematica simplify. *)

(* Express h^2 in terms of rho *)
rho = 2 Cos[Pi/7];

(* Target: h^2 = (2/21)(1 - rho + rho^2) *)
targetHsq = (2/21)(1 - rho + rho^2);

Print["\nStep 4: Compare with Q(rho) formula"];
Print["  Target: h^2 = (2/21)(1 - rho + rho^2)"];
Print["  Target (numeric) = ", N[targetHsq, 30]];

diff = FullSimplify[hSquared - targetHsq];
Print["\n  h^2(derived) - h^2(target) = ", diff];
Print["  Simplified difference = ", FullSimplify[diff, 
  Assumptions -> Element[Pi, Reals]]];

(* Try harder with TrigReduce *)
diff2 = FullSimplify[TrigReduce[hSquared - targetHsq]];
Print["  TrigReduce difference = ", diff2];

(* Numerical check *)
Print["  Numerical difference = ", N[hSquared - targetHsq, 50]];

(* ---- Step 5: Also verify b^2 and h/b ---- *)
Print["\nStep 5: Derived b^2 and h/b"];
bSquared = 8/3 - 7 hSquared;
bSqSimp = FullSimplify[bSquared];
Print["  b^2 = ", N[bSquared, 30]];

targetBsq = 2 + (2/3) rho - (2/3) rho^2;
Print["  Target b^2 = (2 + 2rho/3 - 2rho^2/3) = ", N[targetBsq, 30]];
Print["  b^2 difference (numeric) = ", N[bSquared - targetBsq, 50]];

hOverB = Sqrt[hSquared] / Sqrt[bSquared];
targetHoverB = (-1 + 6 rho - 2 rho^2)/7;
Print["  h/b (numeric) = ", N[hOverB, 30]];
Print["  Target h/b = (-1+6rho-2rho^2)/7 = ", N[targetHoverB, 30]];
Print["  h/b difference (numeric) = ", N[hOverB - targetHoverB, 50]];

(* ---- Step 6: Compute the optimal shadow volume ---- *)
Print["\nStep 6: Optimal shadow volume"];
r = 4/Sqrt[21];
Vzono = r^2 * (Sqrt[hSquared] * S1exact + Sqrt[bSquared] * S2exact);
Vshadow = Vzono / (8/3)^(3/2);
Print["  V_shadow = ", N[Vshadow, 30]];
Print["  Expected: 7.0270375786"];
