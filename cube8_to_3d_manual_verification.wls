rho = 2 Cos[Pi/7];
sigma = Sin[3 Pi/7]/Sin[Pi/7];
northpole = {0, 0, 6 + 9 rho + 4 sigma};
(* internal representation *)
n1 = {2 sigma, -2 rho - 4 sigma, 5 + 7 rho + 3 sigma};
n2 = {4 + 2 rho + 2 sigma, -2 -2 sigma, 5 + 7 rho + 3 sigma};
n3 = {-2 +2 rho + 4 sigma, 2sigma, 5 + 7 rho + 3 sigma};
n4 = {2 + 2 rho - 2 sigma, 4 + 2 rho + 2 sigma, 5 + 7 rho + 3 sigma};
n5 = {-2 - 4 rho, -2 + 2 rho + 4 sigma, 5 + 7 rho + 3 sigma};
n6 = {-2rho -4 sigma, 2+2rho-2sigma, 5 + 7 rho + 3 sigma};
n7 = {-2-2sigma, -2-4rho, 5+ 7 rho + 3 sigma};

(* optimized constants *)
r = 4/Sqrt[21];
h = Sqrt[(2/21)(1 − rho + rho^2)];
b = Sqrt[8/3 - 7 h^2];

shear[internalVector_]:= {internalVector[[1]] + internalVector[[2]]/sigma/2, internalVector[[2]] Sin[3 Pi/7], internalVector[[3]]};

euclideanNodes = Map[ shear, {n1, n2, n3, n4, n5, n6, n7}];

(* verify they are the same *)
Map[Norm, euclideanNodes]//N

(* 2 h *)
(northpole - euclideanNodes)[[3]]

(* distance between two neighbors *)
(* 2 r Sin[Pi/7]*2 *)
EuclideanDistance[euclideanNodes[[1]], euclideanNodes[[2]]]


vZomeRatio = EuclideanDistance[euclideanNodes[[1]], euclideanNodes[[2]]]/(northpole - euclideanNodes[[1]])[[3]];
optimizationRatio = 2 r Sin[Pi/7] * 2/(2 h);

(* Expect them to be the same *)
{vZomeRatio, optimizationRatio}//N

(* 0 *)
FullSimplify[ vZomeRatio - optimizationRatio ]
