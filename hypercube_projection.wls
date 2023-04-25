#!/usr/bin/env wolframscript

(* lowerDimension = 3; *)
lowerDimension = 4;

higherDimension = 2 * lowerDimension;

exportToPov = True;
(* exportToPov = False; *)

(* computeSymbolicMax = True; *)
computeSymbolicMax = False;

(* topView = True; *)
topView = False;

optimal84 = True;
(* optimal84 = False; *)

(* shrinkFactor = 0; *)
shrinkFactor = 0.4;

frameCount = 16;
(* frameCount = 100; *)

lighting = {{"Point", White, {50, -50, 50}}};

imageSize = {4, 3} * 720 / 3;

imageFolder = "results";

shape = "maxVolume_" <> IntegerString[higherDimension];
shapeWithParams = shape <> If[optimal84, "", "_symmetric"] <> "_shrink_"<> IntegerString[ Round[shrinkFactor * 100] ];

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

epsilon = 0.0000001;
randomDirection = {100, Pi, 1, Exp[-4]}//N;
cubeVertices[dimension_] := 2 * Map[IntegerDigits[#, 2][[2;;dimension+1]]&, Table[k, {k, 2^dimension, 2 * 2^dimension - 1}]] - 1;
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];
getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameIds[ids1_, ids2_] := Sort[ids1] == Sort[ids2];

doubleRotationMatrix[theta_] := ArrayFlatten[{{RotationMatrix[-theta], 0}, {0, RotationMatrix[-theta]}}];
(* doubleRotationMatrix[theta_] := ArrayFlatten[{{RotationMatrix[-theta], 0}, {0, RotationMatrix[theta]}}]; *)
(* doubleRotationMatrix[theta_] := ArrayFlatten[{{RotationMatrix[-theta], 0}, {0, IdentityMatrix[2]}}]; *)
projectionMatrix = IdentityMatrix[4][[Range[3]]];

getProjection4[x_, y_] := Module[
	{sx, sy, cx, cy},
	sx = Sin[x/2];
	sy = Sin[y/2];
	cx = Cos[x/2];
	cy = Cos[y/2];
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
(* Map[Print, Transpose[getProjection4[x, y]] ]; *)

phi = (1 + Sqrt[5])/2;
projection3 = {
	{1, phi, 0, -1, phi, 0}, 
	{phi, 0, 1, phi, 0, -1}, 
	{0, 1, phi, 0, -1, phi}
};

orth[target_, ref_] := target - ref * Dot[target, ref]/Dot[ref, ref];
getVolume[vs_] := Module[
	{a, b, c, areaAB, height, bOrth, cOrthA},
	a = vs[[1]];
	b = vs[[2]];
	c = vs[[3]];

	bOrth = orth[b, a];
	areaAB = Norm[a] * Norm[bOrth];
	cOrthA = orth[c, a];
	height = orth[cOrthA, bOrth];

	(* Print[Round[N[{Dot[a, b], Dot[a, bOrth], Dot[cOrthA, a], Dot[height, a], Dot[height, b], Dot[height, bOrth] }], epsilon]]; *)
	areaAB * Norm[height]
];

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

getCellCenter[cid_] := Module[
	{cell, faces4d, faceCenters, cellCenter},
	cell = cells[[cid]];
	faces4d = Map[projectedVertices[[#]]&, cell, {2}];
	faceCenters = Map[Mean, faces4d];
	cellCenter = Mean[faceCenters];
	cellCenter
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
	obj2[x_, y_] := Expand[2 objective[x, y] - 4];
	Print[ obj2[x, y] ];
	(* need to maximize this:
		Sin[x] + 4*Cos[x] + 3*Sin[y] + 4*Cos[y] + 3*Cos[x]*Cos[y]
		numerical max:
		11.6893756408159
	 *)

	Print[D[obj2[x, y], x]];
	Print[D[obj2[x, y], y]];

	(* solution = Solve[ Cos[x] - 4*Sin[x] - 3*Cos[y]*Sin[x] == 0 && 3*Cos[y] - 4*Sin[y] - 3*Cos[x]*Sin[y] == 0, {x, y}, Reals, InverseFunctions -> True]; *)
	(* Print[solution]; *)
	(* still cannot find closed form expression *)

	solution=FindInstance[Cos[x] - 4*Sin[x] - 3*Cos[y]*Sin[x] == 0 && 3*Cos[y] - 4*Sin[y] - 3*Cos[x]*Sin[y] == 0 && x>0 && y>0 && x<1 && y<1, {x, y}]//FullSimplify;

	Print[{x, y}/.solution[[1]]//FullSimplify ];
	Print[{x, y}/.solution[[1]]//N ];
	Print[objective[x, y]/.solution[[1]]//FullSimplify ];
	Print[objective[x, y]/.solution[[1]]//N ];

	Print[{Sin[x/2], Cos[x/2], Sin[y/2], Cos[y/2]}/.solution[[1]]//FullSimplify];
	Print[{Sin[x/2], Cos[x/2], Sin[y/2], Cos[y/2]}/.solution[[1]]//N];

	Exit[];

	maxSolution = Maximize[ {objective[x, y], x>0&&y>0&&x<4/10&&y<6/10}, {x, y} ]//FullSimplify;
	Print[maxSolution];

	Print["maximum"];
	maximum = maxSolution[[1]];
	Print[maximum];
	Print[N[maximum]];
	Print[ToRadicals[maximum]];

	Print["optimal parameters"];
	optimum = {x, y}/.maxSolution[[2]]//FullSimplify;
	Print[optimum];
	Print[N[optimum]];
	Print[ToRadicals[optimum]];

	Exit[];
];

(* back to numeric *)
optimalProjection = getProjection4[x, y]/.maxSolutionN[[2]];
symmetricProjection = getProjection4[x, y]/.{x -> Pi/4, y -> Pi/4};

(* projectedVertices = vertices.Transpose[optimalProjection]/2; *)
projectedVertices = vertices.Transpose[If[optimal84, optimalProjection, symmetricProjection]]/2;

(* calculated in python *)
convexHullIndices = {0, 1, 2, 3, 4, 5, 8, 10, 12, 16, 17, 18, 19, 24, 25, 26, 27, 32, 33, 34, 35, 36, 37, 38, 39, 49, 50, 51, 53, 55, 58, 59, 64, 65, 68, 69, 72, 73, 76, 77, 80, 81, 83, 85, 88, 89, 90, 91, 92, 93, 95, 97, 100, 101, 103, 109, 113, 115, 117, 119, 121, 123, 125, 127, 128, 130, 132, 134, 136, 138, 140, 142, 146, 152, 154, 155, 158, 160, 162, 163, 164, 165, 166, 167, 170, 172, 174, 175, 178, 179, 182, 183, 186, 187, 190, 191, 196, 197, 200, 202, 204, 205, 206, 216, 217, 218, 219, 220, 221, 222, 223, 228, 229, 230, 231, 236, 237, 238, 239, 243, 245, 247, 250, 251, 252, 253, 254, 255} + 1;

vertexIdsInCells = {{0, 4, 8, 12, 64, 68, 72, 76}, {0, 4, 8, 12, 128, 132, 136, 140}, {4, 12, 68, 76, 132, 140, 196, 204}, {8, 12, 72, 76, 136, 140, 200, 204}, {4, 36, 68, 100, 132, 164, 196, 228}, {8, 24, 72, 88, 136, 152, 200, 216}, {64, 65, 68, 69, 72, 73, 76, 77}, {128, 130, 132, 134, 136, 138, 140, 142}, {0, 1, 4, 5, 64, 65, 68, 69}, {0, 2, 8, 10, 128, 130, 136, 138}, {72, 76, 88, 92, 200, 204, 216, 220}, {132, 140, 164, 172, 196, 204, 228, 236}, {0, 4, 32, 36, 128, 132, 160, 164}, {136, 138, 140, 142, 200, 202, 204, 206}, {0, 8, 16, 24, 64, 72, 80, 88}, {68, 69, 76, 77, 196, 197, 204, 205}, {4, 5, 36, 37, 68, 69, 100, 101}, {8, 10, 24, 26, 136, 138, 152, 154}, {72, 73, 76, 77, 88, 89, 92, 93}, {132, 134, 140, 142, 164, 166, 172, 174}, {64, 65, 72, 73, 80, 81, 88, 89}, {128, 130, 132, 134, 160, 162, 164, 166}, {136, 138, 152, 154, 200, 202, 216, 218}, {68, 69, 100, 101, 196, 197, 228, 229}, {0, 1, 16, 17, 64, 65, 80, 81}, {0, 2, 32, 34, 128, 130, 160, 162}, {196, 197, 204, 205, 228, 229, 236, 237}, {200, 202, 204, 206, 216, 218, 220, 222}, {0, 1, 4, 5, 32, 33, 36, 37}, {0, 2, 8, 10, 16, 18, 24, 26}, {76, 77, 92, 93, 204, 205, 220, 221}, {140, 142, 172, 174, 204, 206, 236, 238}, {65, 69, 73, 77, 81, 85, 89, 93}, {130, 134, 138, 142, 162, 166, 170, 174}, {36, 37, 100, 101, 164, 165, 228, 229}, {24, 26, 88, 90, 152, 154, 216, 218}, {16, 17, 24, 25, 80, 81, 88, 89}, {69, 77, 101, 109, 197, 205, 229, 237}, {32, 34, 36, 38, 160, 162, 164, 166}, {138, 142, 154, 158, 202, 206, 218, 222}, {1, 5, 33, 37, 65, 69, 97, 101}, {88, 89, 92, 93, 216, 217, 220, 221}, {164, 166, 172, 174, 228, 230, 236, 238}, {2, 10, 18, 26, 130, 138, 146, 154}, {0, 1, 2, 3, 16, 17, 18, 19}, {0, 1, 2, 3, 32, 33, 34, 35}, {204, 205, 220, 221, 236, 237, 252, 253}, {204, 206, 220, 222, 236, 238, 252, 254}, {16, 17, 18, 19, 24, 25, 26, 27}, {32, 33, 34, 35, 36, 37, 38, 39}, {77, 93, 109, 125, 205, 221, 237, 253}, {142, 158, 174, 190, 206, 222, 238, 254}, {24, 25, 26, 27, 88, 89, 90, 91}, {36, 37, 38, 39, 164, 165, 166, 167}, {69, 77, 85, 93, 101, 109, 117, 125}, {138, 142, 154, 158, 170, 174, 186, 190}, {65, 69, 81, 85, 97, 101, 113, 117}, {130, 138, 146, 154, 162, 170, 178, 186}, {88, 89, 90, 91, 216, 217, 218, 219}, {164, 165, 166, 167, 228, 229, 230, 231}, {2, 18, 34, 50, 130, 146, 162, 178}, {1, 17, 33, 49, 65, 81, 97, 113}, {216, 217, 218, 219, 220, 221, 222, 223}, {228, 229, 230, 231, 236, 237, 238, 239}, {1, 3, 17, 19, 33, 35, 49, 51}, {2, 3, 18, 19, 34, 35, 50, 51}, {220, 221, 222, 223, 252, 253, 254, 255}, {236, 237, 238, 239, 252, 253, 254, 255}, {17, 19, 25, 27, 81, 83, 89, 91}, {34, 35, 38, 39, 162, 163, 166, 167}, {101, 109, 117, 125, 229, 237, 245, 253}, {154, 158, 186, 190, 218, 222, 250, 254}, {33, 37, 49, 53, 97, 101, 113, 117}, {18, 26, 50, 58, 146, 154, 178, 186}, {89, 91, 93, 95, 217, 219, 221, 223}, {166, 167, 174, 175, 230, 231, 238, 239}, {37, 39, 101, 103, 165, 167, 229, 231}, {26, 27, 90, 91, 154, 155, 218, 219}, {81, 85, 89, 93, 113, 117, 121, 125}, {162, 166, 170, 174, 178, 182, 186, 190}, {17, 19, 49, 51, 81, 83, 113, 115}, {34, 35, 50, 51, 162, 163, 178, 179}, {218, 219, 222, 223, 250, 251, 254, 255}, {229, 231, 237, 239, 245, 247, 253, 255}, {18, 19, 26, 27, 50, 51, 58, 59}, {33, 35, 37, 39, 49, 51, 53, 55}, {93, 95, 125, 127, 221, 223, 253, 255}, {174, 175, 190, 191, 238, 239, 254, 255}, {37, 39, 53, 55, 101, 103, 117, 119}, {26, 27, 58, 59, 154, 155, 186, 187}, {89, 91, 93, 95, 121, 123, 125, 127}, {166, 167, 174, 175, 182, 183, 190, 191}, {81, 83, 89, 91, 113, 115, 121, 123}, {162, 163, 166, 167, 178, 179, 182, 183}, {101, 103, 117, 119, 229, 231, 245, 247}, {154, 155, 186, 187, 218, 219, 250, 251}, {49, 51, 53, 55, 113, 115, 117, 119}, {50, 51, 58, 59, 178, 179, 186, 187}, {91, 95, 123, 127, 219, 223, 251, 255}, {167, 175, 183, 191, 231, 239, 247, 255}, {19, 27, 51, 59, 83, 91, 115, 123}, {35, 39, 51, 55, 163, 167, 179, 183}, {117, 119, 125, 127, 245, 247, 253, 255}, {186, 187, 190, 191, 250, 251, 254, 255}, {178, 179, 182, 183, 186, 187, 190, 191}, {39, 55, 103, 119, 167, 183, 231, 247}, {27, 59, 91, 123, 155, 187, 219, 251}, {113, 115, 117, 119, 121, 123, 125, 127}, {51, 55, 115, 119, 179, 183, 243, 247}, {51, 59, 115, 123, 179, 187, 243, 251}, {115, 119, 123, 127, 243, 247, 251, 255}, {179, 183, 187, 191, 243, 247, 251, 255}};
vertexIdsInCells = Map[# + 1 &, vertexIdsInCells, {2}];

convexHullVertices = projectedVertices[[convexHullIndices]];
convexHullVertices = SortBy[convexHullVertices, Last];

cells = {};
cellVolumes = {};

For[cid = 1, cid <= Length[vertexIdsInCells], cid++,
	vertexIdsInCell = vertexIdsInCells[[cid]];
	verticesInCell = projectedVertices[[vertexIdsInCell]];

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
	(* Print[ {cid, tallyTally} ]; *)
	(* {{4, 3}, {2, 6}, {1, 4}}. affine cube, Parallelepiped *)

	cubeEdges = Map[First, Select[tally, #[[2]]==4 &]];
	cubeFaceDiagonals = Map[First, Select[tally, #[[2]]==2 &]];
	(* Print[Map[Round[Norm[#]*Sqrt[2], epsilon]&, cubeEdges]]; *)

	volume = getVolume[cubeEdges];
	dots = {
		Dot[cubeEdges[[1]], cubeEdges[[2]]],
		Dot[cubeEdges[[2]], cubeEdges[[3]]],
		Dot[cubeEdges[[3]], cubeEdges[[1]]]
	};
	dots = Sort[Map[Abs, dots]//N];
	(* AppendTo[cellVolumes, Join[{volume}, dots]]; *)
	(* if we include dots, we are getting the same groups as using only volume *)
	AppendTo[cellVolumes, volume];

	edgeCandidates = Subsets[vertexIdsInCell, {2}];
	edges = {};
	For[eid = 1, eid <= Length[edgeCandidates], eid++,
		edge = edgeCandidates[[eid]];
		edgeCoordinates = projectedVertices[[edge]];
		edgeVector = edgeCoordinates[[1]] - edgeCoordinates[[2]];
		If[edgeVector.randomDirection < 0, edgeVector = - edgeVector];
		If[Length[Intersection[{edgeVector}, cubeEdges, SameTest -> ApproxSamePoint ]] > 0,
			AppendTo[edges, edge];
		];
	];
	edgeGraph = Graph[edges];
	cycles = FindCycle[edgeGraph, {4}, All];
	faces = Map[#[[1]]&, cycles, {2}];
	edgesFromFaces = Union[ Flatten[Map[getEdgesFromFace, faces], 1], SameTest -> sameIds];

	AppendTo[cells, faces];
];

cellVolumesTally = Sort[Tally[cellVolumes//N, ApproxSamePoint[#1, #2]&]];
distinctVolumes = Map[First, cellVolumesTally];
(* Map[Print, cellVolumesTally]; *)

cellIdsByVolume = GroupBy[ Range[Length[cellVolumes]], Round[ cellVolumes[[#]], epsilon]& ];
cellIdGroups = Values[cellIdsByVolume];
cellCentersByGroup = Map[ getCellCenter , cellIdGroups, {2}];

(* export to analyze in python *)
Print[ Map[# - 1 &, cellIdGroups, {2}]];
(* 
In each group, there are 16 cell centers. We treat them as vertices (dualize). 
In each group, the convex hull of 16 points have 56 equations (cells)
Each cell has 4 faces ( I guessed the combination)
In total 56*4/2 = 112 faces

Each cell has 6 edges ( I guessed the combination)
56*6/2 = 336 before dedup
72 unique edges. 4.6667 cells around each edge, on average. Not regular.
4-fold edges x 32
5-fold edges x 32
6-fold edges x 8

Each cell has 4 points (tetrahedron?)
56*4/16 = 14 cells per point


 *)
Exit[];

For[gid =1, gid<=Length[cellCentersByGroup], gid++,
	Print[gid];
	cellCentersInGroup = cellCentersByGroup[[gid]];
	diffs = Map[#[[1]]-#[[2]]&, Subsets[cellCentersInGroup, {2}]];
	diffs = Map[If[ #.randomDirection < 0, -#, # ]&, diffs];

	tally = Tally[diffs, ApproxSamePoint[#1, #2]&];
	(* tally = Tally[Map[Norm, diffs], ApproxSamePoint[#1, #2]&]; *)
	tallyCounts = Map[Last, tally];
	tallyTally = Reverse[Sort[Tally[tallyCounts]]];

	Print[tally];
	Print[tallyTally];
	(* Exit[]; *)
];


Exit[];

Print[cellVolumesTally//Length];
(* Exit[]; *)
colorMap = <||>;
colors = {Red, Blue, Green, Orange, Yellow, Purple, Brown};
For[colorId = 1, colorId <= Length[distinctVolumes], colorId++,
	(* colorMap[Round[ distinctVolumes[[colorId]], epsilon]] = Hue[ colorId/Length[distinctVolumes] ]; *)
	colorMap[Round[ distinctVolumes[[colorId]], epsilon]] = colors[[colorId]];
];

(* Print[colorMap]; *)

(* 
when optimal = false (symmetric)
4 different volumes
{{0.0956708580912724, 16}, {0.23096988312782166, 16}, {0.24999999999999994, 16}, {0.30618621784789724, 64}}

when optimal = true
7 different volumes

{{0.20855353945134233, 16}, {0.24299911162995522, 16}, {0.25266611333566996, 16}, {0.26883897819123576, 16}, {0.2705365471094385, 16}, {0.2725194114324201, 16}, {0.3374691674463592, 16}}
 *)

exportFrame[spinAngle_, frameIndex_, outputFolderName_] := Module[
	{},
	(* If[frameIndex < 77, Return[]]; *)
	(* If[Mod[frameIndex, 2] == 1, Return[]]; *)
	(* Return[]; *)

	elements = {};

	For[cid = 1, cid <= Length[vertexIdsInCells], cid++,
		cell = cells[[cid]];
		faces4d = Map[projectedVertices[[#]]&, cell, {2}];
		faceCenters = Map[Mean, faces4d];
		cellCenter = Mean[faceCenters];

		shrunkFaces4d = Map[ shrinkFactor * cellCenter + (1 - shrinkFactor) *# &, faces4d, {2} ];

		(* Print[cellCenter]; *)
		(* Exit[]; *)
		(* Print[faces4d]; *)
		faces3d = shrunkFaces4d.doubleRotationMatrix[spinAngle].Transpose[projectionMatrix];
		edges3d = Apply[Join, Map[getEdgesFromFace, faces3d]];
		color = colorMap[ Round[cellVolumes[[cid]], epsilon] ]; 

		thickness = 0.003;
		AppendTo[elements, {Opacity[1], FaceForm[color], Map[ Polygon, faces3d ]}];

		(* AppendTo[elements, {Opacity[1], Map[ Polygon, faces3d ]}]; *)

		AppendTo[elements, {Opacity[1], Map[ Cylinder[#, thickness]&, edges3d ]}];
	];

	(* graph = Plot3D[getProjectedVolume[getProjection4[x, y]], {x, 0, Pi/4}, {y, 0, Pi/4}]; *)
	(* graph = ContourPlot[Log[(target4 - getProjectedVolume[getProjection4[x/180*Pi, y/180*Pi]])], {x, 0, 45}, {y, 0, 45}, Contours -> 100]; *)

	(* graph = ContourPlot[Exp[objective[x/180*Pi, y/180*Pi]], {x, 0, 90}, {y, 0, 90}, Contours -> 100]; *)

	graph = Graphics3D[
			elements,
			ViewPoint -> If[topView, {0, 0, 1000}, {1000, -1000, 500}],
			ViewVertical -> If[topView, {0, 1, 0}, {0, 0, 1}],
			Lighting -> lighting,
			Boxed -> False
		];

	(* outputFolderName = "hypercube"; *)
	outputFolder = FileNameJoin[{imageFolder, outputFolderName}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	outputFileName = "spinning_frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
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
];

exportAnimation[] := Module[
	{
		outputFolder,
		gifFileName,
		gifFolder,
		maxExplosionFactor,
		explodeStep,
		explosionFactors,
		frameIndex,
		explosionFactor,
		gifCommand,
		edges
	},

	outputFolderName = shapeWithParams <> "_spinning_frames";
	outputFolder = FileNameJoin[{imageFolder, outputFolderName}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	gifFileName = shapeWithParams <> "_spinning";
	gifFolder = FileNameJoin[{imageFolder, shape}];
	If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];

	For[frameIndex = 1, frameIndex <= frameCount, frameIndex++,	
		spinAngle = frameIndex/frameCount * Pi;
		Print[DateString["ISODateTime"] <> " frameIndex/frameCount: " <> IntegerString[frameIndex] <> " / " <> IntegerString[frameCount] <> " spin angle/Pi: " <> ToString[spinAngle/Pi//N] ];

		exportFrame[spinAngle, frameIndex, outputFolderName];
	];

	(* Exit[]; *)

	delay = If[frameCount >= 40, 3, Floor[300/frameCount]];
	(* delay = 2; *)

	gifCommand = "magick convert -quiet +repage -alpha remove -loop 0 -delay " <> IntegerString[delay] <> " -dispose 2 " <> outputFolder <> "\\spinning_frame*.png " <> gifFolder <> "\\" <> gifFileName <> ".gif";

	Print[gifCommand];
	Run[gifCommand];
	Print["Saved the animation to " <> gifFolder <> "\\" <> gifFileName <> ".gif"];
	
  mp4Command = "ffmpeg -i " <> gifFolder <> "\\" <> gifFileName <> ".gif" 
      <> " -movflags faststart -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" "
      <> gifFolder <> "\\" <> gifFileName <> ".mp4";

  mp4FileName = gifFolder <> "\\" <> gifFileName <> ".mp4";
  If[FileExistsQ[mp4FileName], DeleteFile[mp4FileName]];
  Print[mp4Command];
  Run[mp4Command];
  Print["Saved the video to " <> gifFolder <> "\\" <> gifFileName <> ".mp4"];
];

exportAnimation[];


(* For[id=1, id<=10, id++,
	exportFrame[id*Pi/2, id, "test"];
];
 *)
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
{
2*ArcTan[Root[29 - 464*#1 + 932*#1^2 + 208*#1^3 + 774*#1^4 + 176*#1^5 + 308*#1^6 + 80*#1^7 + 5*#1^8 & , 3, 0]], 
ArcCos[Root[4096 + 10752*#1 - 28928*#1^2 - 107040*#1^3 - 78272*#1^4 + 58920*#1^5 + 116712*#1^6 + 62640*#1^7 + 11745*#1^8 & , 4, 0]]}

{0.14696147003406254, 0.406571006373532}


 *)