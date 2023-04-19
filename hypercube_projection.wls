#!/usr/bin/env wolframscript

(* lowerDimension = 3; *)
lowerDimension = 4;

higherDimension = 2 * lowerDimension;

(* exportToPov = True; *)
exportToPov = False;

(* computeSymbolicMax = True; *)
computeSymbolicMax = False;

imageSize = {4, 3} * 720 / 3;

imageFolder = "results";


Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];
cubeVertices[dimension_] := 2 * Map[IntegerDigits[#, 2][[2;;dimension+1]]&, Table[k, {k, 2^dimension, 2 * 2^dimension - 1}]] - 1;
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];

epsilon = 0.00000000001;
randomDirection = {100, Pi, 1, Exp[-4]}//N;

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
	{norm, normalizedProjection, subsets},
	(* this is 4 * identity matrix *)
	(* Print[projection.Transpose[projection]//Simplify]; *)
	norm = If[lowerDimension == 4, 2, Norm[projection[[1]]]//Simplify];
	normalizedProjection = Transpose[projection/norm];
	subsets = Subsets[normalizedProjection, {lowerDimension}];
	Total[Map[Abs[Det[#]]&, subsets]]
];

getSigns[projection_] := Module[
	{norm, normalizedProjection, subsets},
	norm = If[lowerDimension == 4, 2, Norm[projection[[1]]]//Simplify];
	normalizedProjection = Transpose[projection/norm];
	subsets = Subsets[normalizedProjection, {lowerDimension}];
	Map[Sign[Det[#]]&, subsets]
];

getProjectedVolumeBySigns[projection_, signs_] := Module[
	{norm, normalizedProjection, subsets},
	(* this is 4 * identity matrix *)
	(* Print[projection.Transpose[projection]//Simplify]; *)
	norm = If[lowerDimension == 4, 2, Norm[projection[[1]]]//Simplify];
	normalizedProjection = Transpose[projection/norm];
	subsets = Subsets[normalizedProjection, {lowerDimension}];
	Total[Table[ signs[[k]] * Det[subsets[[k]] ] , {k, Length[subsets]}]]
];


optimalX = 0.07348073502590648;
optimalY = 0.20328550320237618;
projection4 = getProjection4[optimalX, optimalY];
(* Print[projection4]; *)

projection = If[lowerDimension == 3, projection3, projection4];

vertices = Map[# + Table[0, {k, higherDimension}]&, cubeVertices[higherDimension]]/2;
(* Print[vertices]; *)

(* Print[projection3//N]; *)

(* projectedVertices = vertices.Transpose[Orthogonalize[projection]]; *)
(* Print[projectedVertices//Length]; *)

maxSolutionN = NMaximize[ {getProjectedVolume[getProjection4[x, y]], x>0 && y>0}, {x, y} ];
Print[maxSolutionN];
optimumN = {x, y}/.maxSolutionN[[2]];
Print["optimum"];
Print[optimumN];
Print[optimumN*180/Pi];

If[computeSymbolicMax,
	signs = getSigns[Apply[getProjection4, optimumN]];
	objective[x_, y_] := getProjectedVolumeBySigns[getProjection4[x, y], signs]//FullSimplify;
	maxSolution = Maximize[ {objective[x, y], x>0&&y>0&&x<2/10&&y<3/10}, {x, y} ]//FullSimplify;
	Print[maxSolution];

	Print["maximum"];
	maximum = maxSolution[[1]];
	Print[maximum];
	Print[N[maximum]];

	Print["optimal parameters"];
	optimum = {x, y}/.maxSolution[[2]]//FullSimplify;
	Print[optimum];
	Print[N[optimum]];
];

(* back to numeric *)
optimalProjection = getProjection4[x, y]/.maxSolutionN[[2]];

projectedVertices = vertices.Transpose[optimalProjection]/2;

(* calculated in python *)
convexHullIndices = {0, 1, 2, 3, 4, 5, 8, 10, 12, 16, 17, 18, 19, 24, 25, 26, 27, 32, 33, 34, 35, 36, 37, 38, 39, 49, 50, 51, 53, 55, 58, 59, 64, 65, 68, 69, 72, 73, 76, 77, 80, 81, 83, 85, 88, 89, 90, 91, 92, 93, 95, 97, 100, 101, 103, 109, 113, 115, 117, 119, 121, 123, 125, 127, 128, 130, 132, 134, 136, 138, 140, 142, 146, 152, 154, 155, 158, 160, 162, 163, 164, 165, 166, 167, 170, 172, 174, 175, 178, 179, 182, 183, 186, 187, 190, 191, 196, 197, 200, 202, 204, 205, 206, 216, 217, 218, 219, 220, 221, 222, 223, 228, 229, 230, 231, 236, 237, 238, 239, 243, 245, 247, 250, 251, 252, 253, 254, 255} + 1;

vertexIdsInCells = {{0, 4, 8, 12, 64, 68, 72, 76}, {0, 4, 8, 12, 128, 132, 136, 140}, {4, 12, 68, 76, 132, 140, 196, 204}, {8, 12, 72, 76, 136, 140, 200, 204}, {4, 36, 68, 100, 132, 164, 196, 228}, {8, 24, 72, 88, 136, 152, 200, 216}, {64, 65, 68, 69, 72, 73, 76, 77}, {128, 130, 132, 134, 136, 138, 140, 142}, {0, 1, 4, 5, 64, 65, 68, 69}, {0, 2, 8, 10, 128, 130, 136, 138}, {72, 76, 88, 92, 200, 204, 216, 220}, {132, 140, 164, 172, 196, 204, 228, 236}, {0, 4, 32, 36, 128, 132, 160, 164}, {136, 138, 140, 142, 200, 202, 204, 206}, {0, 8, 16, 24, 64, 72, 80, 88}, {68, 69, 76, 77, 196, 197, 204, 205}, {4, 5, 36, 37, 68, 69, 100, 101}, {8, 10, 24, 26, 136, 138, 152, 154}, {72, 73, 76, 77, 88, 89, 92, 93}, {132, 134, 140, 142, 164, 166, 172, 174}, {64, 65, 72, 73, 80, 81, 88, 89}, {128, 130, 132, 134, 160, 162, 164, 166}, {136, 138, 152, 154, 200, 202, 216, 218}, {68, 69, 100, 101, 196, 197, 228, 229}, {0, 1, 16, 17, 64, 65, 80, 81}, {0, 2, 32, 34, 128, 130, 160, 162}, {196, 197, 204, 205, 228, 229, 236, 237}, {200, 202, 204, 206, 216, 218, 220, 222}, {0, 1, 4, 5, 32, 33, 36, 37}, {0, 2, 8, 10, 16, 18, 24, 26}, {76, 77, 92, 93, 204, 205, 220, 221}, {140, 142, 172, 174, 204, 206, 236, 238}, {65, 69, 73, 77, 81, 85, 89, 93}, {130, 134, 138, 142, 162, 166, 170, 174}, {36, 37, 100, 101, 164, 165, 228, 229}, {24, 26, 88, 90, 152, 154, 216, 218}, {16, 17, 24, 25, 80, 81, 88, 89}, {69, 77, 101, 109, 197, 205, 229, 237}, {32, 34, 36, 38, 160, 162, 164, 166}, {138, 142, 154, 158, 202, 206, 218, 222}, {1, 5, 33, 37, 65, 69, 97, 101}, {88, 89, 92, 93, 216, 217, 220, 221}, {164, 166, 172, 174, 228, 230, 236, 238}, {2, 10, 18, 26, 130, 138, 146, 154}, {0, 1, 2, 3, 16, 17, 18, 19}, {0, 1, 2, 3, 32, 33, 34, 35}, {204, 205, 220, 221, 236, 237, 252, 253}, {204, 206, 220, 222, 236, 238, 252, 254}, {16, 17, 18, 19, 24, 25, 26, 27}, {32, 33, 34, 35, 36, 37, 38, 39}, {77, 93, 109, 125, 205, 221, 237, 253}, {142, 158, 174, 190, 206, 222, 238, 254}, {24, 25, 26, 27, 88, 89, 90, 91}, {36, 37, 38, 39, 164, 165, 166, 167}, {69, 77, 85, 93, 101, 109, 117, 125}, {138, 142, 154, 158, 170, 174, 186, 190}, {65, 69, 81, 85, 97, 101, 113, 117}, {130, 138, 146, 154, 162, 170, 178, 186}, {88, 89, 90, 91, 216, 217, 218, 219}, {164, 165, 166, 167, 228, 229, 230, 231}, {2, 18, 34, 50, 130, 146, 162, 178}, {1, 17, 33, 49, 65, 81, 97, 113}, {216, 217, 218, 219, 220, 221, 222, 223}, {228, 229, 230, 231, 236, 237, 238, 239}, {1, 3, 17, 19, 33, 35, 49, 51}, {2, 3, 18, 19, 34, 35, 50, 51}, {220, 221, 222, 223, 252, 253, 254, 255}, {236, 237, 238, 239, 252, 253, 254, 255}, {17, 19, 25, 27, 81, 83, 89, 91}, {34, 35, 38, 39, 162, 163, 166, 167}, {101, 109, 117, 125, 229, 237, 245, 253}, {154, 158, 186, 190, 218, 222, 250, 254}, {33, 37, 49, 53, 97, 101, 113, 117}, {18, 26, 50, 58, 146, 154, 178, 186}, {89, 91, 93, 95, 217, 219, 221, 223}, {166, 167, 174, 175, 230, 231, 238, 239}, {37, 39, 101, 103, 165, 167, 229, 231}, {26, 27, 90, 91, 154, 155, 218, 219}, {81, 85, 89, 93, 113, 117, 121, 125}, {162, 166, 170, 174, 178, 182, 186, 190}, {17, 19, 49, 51, 81, 83, 113, 115}, {34, 35, 50, 51, 162, 163, 178, 179}, {218, 219, 222, 223, 250, 251, 254, 255}, {229, 231, 237, 239, 245, 247, 253, 255}, {18, 19, 26, 27, 50, 51, 58, 59}, {33, 35, 37, 39, 49, 51, 53, 55}, {93, 95, 125, 127, 221, 223, 253, 255}, {174, 175, 190, 191, 238, 239, 254, 255}, {37, 39, 53, 55, 101, 103, 117, 119}, {26, 27, 58, 59, 154, 155, 186, 187}, {89, 91, 93, 95, 121, 123, 125, 127}, {166, 167, 174, 175, 182, 183, 190, 191}, {81, 83, 89, 91, 113, 115, 121, 123}, {162, 163, 166, 167, 178, 179, 182, 183}, {101, 103, 117, 119, 229, 231, 245, 247}, {154, 155, 186, 187, 218, 219, 250, 251}, {49, 51, 53, 55, 113, 115, 117, 119}, {50, 51, 58, 59, 178, 179, 186, 187}, {91, 95, 123, 127, 219, 223, 251, 255}, {167, 175, 183, 191, 231, 239, 247, 255}, {19, 27, 51, 59, 83, 91, 115, 123}, {35, 39, 51, 55, 163, 167, 179, 183}, {117, 119, 125, 127, 245, 247, 253, 255}, {186, 187, 190, 191, 250, 251, 254, 255}, {178, 179, 182, 183, 186, 187, 190, 191}, {39, 55, 103, 119, 167, 183, 231, 247}, {27, 59, 91, 123, 155, 187, 219, 251}, {113, 115, 117, 119, 121, 123, 125, 127}, {51, 55, 115, 119, 179, 183, 243, 247}, {51, 59, 115, 123, 179, 187, 243, 251}, {115, 119, 123, 127, 243, 247, 251, 255}, {179, 183, 187, 191, 243, 247, 251, 255}};
vertexIdsInCells = Map[# + 1 &, vertexIdsInCells, {2}];

convexHullVertices = projectedVertices[[convexHullIndices]];
convexHullVertices = SortBy[convexHullVertices, Last];

verticesInCells = Map[
	( projectedVertices[[#]] )&
	, vertexIdsInCells
];
Print[Length[verticesInCells]];

For[cid = 1, cid <= Length[verticesInCells], cid++,
	verticesInCell = verticesInCells[[cid]];

	(* cellCenter = Mean[verticesInCell]; *)
	(* relative = Map[(#-cellCenter)&, verticesInCell]; *)
	(* negativeRelative = Map[(cellCenter-#)&, verticesInCell]; *)
	(* Print[Length[Union[Join[relative, negativeRelative], SameTest -> ApproxSamePoint]]]; *)

	(* Print[verticesInCell//Length]; *)
	(* Print[verticesInCell]; *)

	diffs = Map[#[[1]]-#[[2]]&, Subsets[verticesInCell, {2}]];
	diffs = Map[If[ #.randomDirection < 0, -#, # ]&, diffs];

	tally = Tally[diffs, ApproxSamePoint[#1, #2]&];
	tallyCounts = Map[Last, tally];
	tallyTally = Reverse[Sort[Tally[tallyCounts]]];
	Print[ {cid, tallyTally} ];
	(* {{4, 3}, {2, 6}, {1, 4}}. affine cube, Parallelepiped *)

	(* Exit[]; *)
];
(* Print[tally]; *)


Exit[];

(* Print[vertexIdsInCells]; *)


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

maximum
Root[-6880 + 13420*#1 + 14883*#1^2 + 712*#1^3 - 2597*#1^4 - 1572*#1^5 - 1359*#1^6 - 432*#1^7 + 81*#1^8 & , 4, 0]
7.844687820407956

optimal parameters
{ArcTan[Root[{-32820224 + 63857664*#1 - 40975680*#1^2 + 12708864*#1^3 - 2081552*#1^4 + 143040*#1^5 + 6660*#1^6 - 1728*#1^7 + 81*#1^8 & , -42 + 8*#1 + #1^2 + 16*#2 + 4*#1*#2 - 60*#2^2 + 2*#1^2*#2^2 - 16*#2^3 + 4*#1*#2^3 + 6*#2^4 - 8*#1*#2^4 + #1^2*#2^4 & }, {1, 2}]], ArcTan[Root[{-32820224 + 63857664*#1 - 40975680*#1^2 + 12708864*#1^3 - 2081552*#1^4 + 143040*#1^5 + 6660*#1^6 - 1728*#1^7 + 81*#1^8 & , -42 + 8*#1 + #1^2 + 16*#2 + 4*#1*#2 - 60*#2^2 + 2*#1^2*#2^2 - 16*#2^3 + 4*#1*#2^3 + 6*#2^4 - 8*#1*#2^4 + #1^2*#2^4 & , 11 + #1 + 2*#2 - 3*#2^2 + #1*#2^2 + 6*#3 + 6*#2^2*#3 - 3*#3^2 + #1*#3^2 + 2*#2*#3^2 - 5*#2^2*#3^2 + #1*#2^2*#3^2 & }, {1, 2, 2}]]}
{0.07348073501703127, 0.20328550318676594}

in degrees
{4.210135992503563, 11.647401369689337}


 *)