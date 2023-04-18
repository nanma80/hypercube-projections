#!/usr/bin/env wolframscript

(* lowerDimension = 3; *)
lowerDimension = 4;

higherDimension = 2 * lowerDimension;

(* exportToPov = True; *)
exportToPov = False;

imageSize = {4, 3} * 720 / 3;

imageFolder = "results";

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];
cubeVertices[dimension_] := 2 * Map[IntegerDigits[#, 2][[2;;dimension+1]]&, Table[k, {k, 2^dimension, 2 * 2^dimension - 1}]] - 1;
epsilon = 0.00000000001;

getProjection4[x_, y_] := Module[
	{sx, sy, cx, cy},
	sx = Sin[x];
	sy = Sin[y];
	cx = Cos[x];
	cy = Cos[y];
	{
		{sx, sx, cx, cx, -sx, -sx, cx, cx},
		{cx, cx, sx, sx, cx, cx, -sx, -sx},
		{cy, -cy, sy, -sy, sy, -sy, cy, -cy},
		{sy, -sy, cy, -cy, -cy, cy, -sy, sy}
	}
];
(* 
getProjectedVolume[getProjection4[x, y]] is symmetric wrt swapping x <-> y
*)


Map[Print, Transpose[getProjection4[x, y]] ];

phi = (1 + Sqrt[5])/2;
projection3 = {
	{1, phi, 0, -1, phi, 0}, 
	{phi, 0, 1, phi, 0, -1}, 
	{0, 1, phi, 0, -1, phi}
};

getProjectedVolume[projection_] := Module[
	{},
	(* this is 4 * identity matrix *)
	(* Print[projection.Transpose[projection]//Simplify]; *)
	norm = If[lowerDimension ==4, 2, Norm[projection[[1]]]//Simplify];
	normalizedProjection = Transpose[projection/norm];
	subsets = Subsets[normalizedProjection, {lowerDimension}];
	Total[Map[Abs[Det[#]]&, subsets]]
];

optimalX = 0.71191742;
optimalY = 0.20328551;

projection4 = getProjection4[optimalX, optimalY];
(* Print[projection4]; *)

projection = If[lowerDimension == 3, projection3, projection4];

vertices = Map[# + Table[0, {k, higherDimension}]&, cubeVertices[higherDimension]]/2;
(* Print[vertices]; *)

(* Print[projection3//N]; *)

(* projectedVertices = vertices.Transpose[Orthogonalize[projection]]; *)
(* Print[projectedVertices//Length]; *)

projectedVolume = getProjectedVolume[projection];
Print["preset value"];
Print[projectedVolume];
Print[projectedVolume//N];

maxSolution = NMaximize[ {getProjectedVolume[getProjection4[x, y]], x>0 && y>0}, {x, y} ];
Print[maxSolution];
optimum = {x, y, x+y, x-y}/.maxSolution[[2]];
Print["optimum"];
Print[optimum];
Print[optimum*180/Pi];

For[id = 1, id <= 4, id++,
	Print[{id, Cos[optimum[[id]]], Sin[optimum[[id]]], Tan[optimum[[id]]] }];
];

(* Print["4-vectors"]; *)
(* Map[Print, Transpose[getProjection4[x, y]]/.maxSolution[[2]]]; *)

optimalProjection = getProjection4[x, y]/.maxSolution[[2]];
(* Map[Print, optimalProjection]; *)

projectedVertices = vertices.Transpose[optimalProjection];

(* calculated in python *)
convexHullIndices = {0, 1, 2, 3, 4, 5, 8, 10, 12, 16, 17, 18, 19, 24, 25, 26, 27, 32, 33, 34, 35, 36, 37, 38, 39, 49, 50, 51, 53, 55, 58, 59, 64, 65, 68, 69, 72, 73, 76, 77, 80, 81, 83, 85, 88, 89, 90, 91, 92, 93, 95, 97, 100, 101, 103, 109, 113, 115, 117, 119, 121, 123, 125, 127, 128, 130, 132, 134, 136, 138, 140, 142, 146, 152, 154, 155, 158, 160, 162, 163, 164, 165, 166, 167, 170, 172, 174, 175, 178, 179, 182, 183, 186, 187, 190, 191, 196, 197, 200, 202, 204, 205, 206, 216, 217, 218, 219, 220, 221, 222, 223, 228, 229, 230, 231, 236, 237, 238, 239, 243, 245, 247, 250, 251, 252, 253, 254, 255} + 1;

convexHullVertices = projectedVertices[[convexHullIndices]];
convexHullVertices = SortBy[convexHullVertices, Last];

Print["convex hull vertices"];
Map[Print, convexHullVertices];

Print[convexHullVertices//Length];

Exit[];

(* maxSolution2 = NMaximize[ getProjectedVolume[getProjection4[0, y]], {y} ]; *)
(* Print[maxSolution2]; *)


(* calculate volume using ConvexHullRegion, which only works for 3D. No longer needed *)
(* chN = ConvexHullRegion[projectedVertices//N];
volumeN = Volume[chN];
Print[volumeN];


ch = ConvexHullRegion[projectedVertices];
volume = Volume[ch];
Print[volume//FullSimplify];
 *)
(* Print[vertices.Transpose[projection3]//N]; *)


target4 = 7.844687820407957;
(* graph = Plot3D[getProjectedVolume[getProjection4[x, y]], {x, 0, Pi/4}, {y, 0, Pi/4}]; *)
(* graph = ContourPlot[Log[(target4 - getProjectedVolume[getProjection4[x/180*Pi, y/180*Pi]])], {x, 0, 45}, {y, 0, 45}, Contours -> 100]; *)

outputFolderName = "hypercube";
outputFolder = FileNameJoin[{imageFolder, outputFolderName}];
If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

outputFileName = "maxVolume_" <> IntegerString[higherDimension] <> ".png";
outputFileName = FileNameJoin[{outputFolder, outputFileName}];

cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

If[exportToPov,
  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
  Run[cropCommand];
  ,
  Export[outputFileName, graph];
];
Print["Saved an image to " <> outputFileName];

(* 

6->3: phi based vertices.
volume: 
4.352501798965643
Sqrt[10 + 4*Sqrt[5]]

8->4: 
target volume: 
7.844687820407957



 *)