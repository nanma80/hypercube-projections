#!/usr/bin/env wolframscript

(* lowerDimension = 3; *)
lowerDimension = 4;

higherDimension = 2 * lowerDimension;

exportToPov = True;
(* exportToPov = False; *)

imageSize = {4, 3} * 720 / 3;

imageFolder = "results";

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];
cubeVertices[dimension_] := 2 * Map[IntegerDigits[#, 2][[2;;dimension+1]]&, Table[k, {k, 2^dimension, 2 * 2^dimension - 1}]] - 1;

getProjection4[x_, y_] := Module[
	{a, b, c, d},
	a = Cos[x]/2;
	b = Sin[x]/2;
	c = Cos[y]/2;
	d = Sin[y]/2;
	{
		{a, a,  a,  a, b, b, b,  b},
		{b, b, -b, -b, a, a, -a, -a},
		{c, -c, d, -d, d, -d, c, -c},
		{d, -d, c, -c, -c, c, -d, d}
	}
];

phi = (1 + Sqrt[5])/2;
projection3 = {
	{1, phi, 0, -1, phi, 0}, 
	{phi, 0, 1, phi, 0, -1}, 
	{0, 1, phi, 0, -1, phi}
};

optimalX = 0.71191742;  
optimalY = 0.20328551;

projection4 = getProjection4[optimalX, optimalY];
Print[projection4];

projection = If[lowerDimension == 3, projection3, projection4];

vertices = Map[# + Table[1, {k, higherDimension}]&, cubeVertices[higherDimension]]/2;

Print[vertices];

(* Print[projection3//N]; *)

projectedVertices = vertices.Transpose[Orthogonalize[projection]];
Print[projectedVertices//Length];

Exit[];

chN = ConvexHullRegion[projectedVertices//N];
volumeN = Volume[chN];
Print[volumeN];


ch = ConvexHullRegion[projectedVertices];
volume = Volume[ch];
Print[volume//FullSimplify];

(* Print[vertices.Transpose[projection3]//N]; *)



Exit[];
graph = Graphics3D[Sphere[{0, 0, 0}, 1]];


outputFolderName = "hypercube";
outputFolder = FileNameJoin[{imageFolder, outputFolderName}];
If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

outputFileName = "cylinder_" <> IntegerString[higherDimension] <> ".png";
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


 *)