(* Symbolic approach: solve q^2 = 4*alpha*(4-alpha)/63 in Q(rho) *)
(* where alpha = 1-rho+rho^2, rho = 2cos(pi/7), rho^3 = rho^2 + 2rho - 1 *)

(* Let q = p0 + p1*rho + p2*rho^2. *)
(* q^2 expanded mod rho^3 = rho^2+2rho-1: *)
(*   const:  p0^2 - 2*p1*p2 - p2^2 *)
(*   rho:    2*p0*p1 + 4*p1*p2 + p2^2 *)  
(*   rho^2:  p1^2 + 2*p0*p2 + 2*p1*p2 + 3*p2^2 *)

(* 4*alpha*(4-alpha)/63 where alpha = 1-rho+rho^2 *)
(* alpha*(4-alpha) = (1-rho+rho^2)(3+rho-rho^2) *)
(* Expand mod rho^3=rho^2+2rho-1: *)
(* = 3+rho-rho^2 -3rho-rho^2+rho^3 +3rho^2+rho^3-rho^4 *)
(* rho^3 = rho^2+2rho-1, rho^4 = 3rho^2+rho-1 *)
(* = 3+rho-rho^2-3rho-rho^2+(rho^2+2rho-1)+3rho^2+(rho^2+2rho-1)-(3rho^2+rho-1) *)
(* = 3+rho-rho^2-3rho-rho^2+rho^2+2rho-1+3rho^2+rho^2+2rho-1-3rho^2-rho+1 *)
(* = (3-1-1+1) + (1-3+2+2-1)rho + (-1-1+1+3+1-3)rho^2 *)
(* = 2 + rho + 0*rho^2 *)
(* So alpha*(4-alpha) = 2+rho *)
(* h^2*b^2 = 4(2+rho)/63 *)

Print["alpha*(4-alpha) = (1-rho+rho^2)(3+rho-rho^2)"];
Print[];

(* Verify symbolically *)
rho = 2 Cos[Pi/7];
alpha = 1 - rho + rho^2;
check = FullSimplify[alpha (4 - alpha) - (2 + rho)];
Print["alpha*(4-alpha) - (2+rho) = ", check];
Print[];

(* So h^2*b^2 = 4(2+rho)/63 *)
Print["h^2*b^2 = 4(2+rho)/63"];
Print["  = 8/63 + (4/63)*rho"];
Print["  numeric = ", N[4(2 + rho)/63, 20]];
Print["  expected = ", N[(2/21)*(1-rho+rho^2)*(8/3 - (2/3)*(1-rho+rho^2)), 20]];
Print[];

(* Now solve: q^2 = 8/63 + (4/63)*rho + 0*rho^2 *)
(* i.e. *)
(*   p0^2 - 2*p1*p2 - p2^2 = 8/63 *)
(*   2*p0*p1 + 4*p1*p2 + p2^2 = 4/63 *)
(*   p1^2 + 2*p0*p2 + 2*p1*p2 + 3*p2^2 = 0 *)

{p0, p1, p2} = {p0, p1, p2};
eqs = {
  p0^2 - 2 p1 p2 - p2^2 == 8/63,
  2 p0 p1 + 4 p1 p2 + p2^2 == 4/63,
  p1^2 + 2 p0 p2 + 2 p1 p2 + 3 p2^2 == 0
};

Print["Solving q^2 = 4(2+rho)/63 for q = p0 + p1*rho + p2*rho^2:"];
sols = Solve[eqs, {p0, p1, p2}];
Print["Number of solutions: ", Length[sols]];
Print[];

Do[
  {a, b, c} = {p0, p1, p2} /. sol;
  val = N[a + b * 2 Cos[Pi/7] + c * (2 Cos[Pi/7])^2, 15];
  Print["Solution: p0=", a, ", p1=", b, ", p2=", c];
  Print["  Value at rho_0: ", val];
  If[val > 0, Print["  *** POSITIVE ***"]];
  Print[];
, {sol, sols}];
