#!/usr/bin/env python

import datetime
import numpy as np
import sympy as sp
from scipy.spatial import ConvexHull
import polytope_builder as pb

dimension = 4


def convex_hull(points):
  ch = ConvexHull(points)
  ch_vertices = ch.vertices
  return [points[index] for index in ch_vertices]


def zonohedrify(all_vertices):
	print "zonohedrifying ", len(all_vertices), " points"
	points = [ np.array(pb.get_origin(dimension)) ]
	vertex_count = -1
	for v in all_vertices:
		vertex_count += 1
		negative_v = [ - elem for elem in v]
		new_points1 = [[sum(x) for x in zip(p, v)] for p in points]
		new_points2 = [[sum(x) for x in zip(p, negative_v)] for p in points]

		points = np.concatenate((new_points1, new_points2), axis = 0)
		if len(points) > 10 * len(all_vertices):
			points = convex_hull(points)
		print datetime.datetime.now(), vertex_count, len(points)

	points = convex_hull(points)
	points.sort(key=lambda x:x[-1] * 1000000 + x[-2] * 10000 + x[-3] * 100 + x[-4])
	return points


def get_runcinated_5_cell_vertices():
	vertices_5_cell = pb.get_5_cell_vertices()
	vertices = []
	for i in xrange(5):
		for j in xrange(5):
			if i == j:
				continue
			vertex = [x[0] - x[1] for x in zip(vertices_5_cell[i], vertices_5_cell[j])]
			vertices.append(vertex)
	return vertices


has_center_symmetry = True

# all_vertices, has_center_symmetry = pb.get_5_cell_vertices(), False
all_vertices = get_runcinated_5_cell_vertices()
# all_vertices = pb.get_orthoplex_vertices(dimension)
# all_vertices = pb.get_cube_vertices(dimension)
# all_vertices = pb.get_24_cell_vertices()
# all_vertices = pb.get_600_cell_vertices()
# all_vertices = pb.get_120_cell_vertices()

# get all edge centers. Note that edge centers of 16-cell are the vertices of 24-cell
# all_vertices = pb.get_edge_centers(all_vertices)

# optimization to remove redundant vertices. 
# Should not include for 5-cell based points, because they don't have center symmetry
random_direction = np.array([0.123, 0.23, 0.76, 0.771])
if has_center_symmetry:
  all_vertices = [v for v in all_vertices if np.inner(v, random_direction) > 0]



points = zonohedrify(all_vertices)
# print all_vertices

scaling = sp.sqrt(1) / 2
points = [p * scaling for p in points]

print "vertex count:", len(points)
print "vertices:"
np.set_printoptions(threshold=np.nan)
print repr(np.array([list(p) for p in points]))
# for p in points:
# 	print p

# investigating zono(tesseract)
# selected_points = [p for p in points if 
# 	1==1
# 	# and (p[2] + p[3] == 4) 
# 	and (p[0] + p[1] + p[2] + p[3] == 6)
# 	]
# print len(selected_points)
# for p in selected_points:
# 	print p

# 5-cell -> 30 vertices, 20 cells: 20 rhombohedra. Dual of runcinated 5-cell. 4 cells around north pole, projected to a rhombic dodecahedron just like the vertex first projection of tesseract. 4 cells around south pole. 12 cells along equator, projected to flat rhombi. Edges crossing equator directly connect (x,y,z,w) to (x,y,z,w+4/sqrt(5))
# runcinated 5-cell -> omnitruncated 5-cell
# 16-cell -> 16 vertices, hypercube
# 8-cell -> 104 vertices, 40 cells: 24 dihedral rhombic dodecahedron + 16 cubes. Can be seen as the rectified 24-cell with eight cubes augmented by pyramids. A rhombic dodecahedron has eight 60-120 rhombi and four squares
# Since 8-cell has 16 vertices form 8 pairs, this zonotope can be viewed as a projection from 8D -> 4D
# 24-cell -> 192 vertices: truncated 24-cell
# 600-cell -> 7616 vertices
# 120-cell didn't finish, too slow
