import numpy as np
import sympy as sp
from scipy.spatial import ConvexHull
import polytope_builder as pb

def convex_hull(points):
	ch = ConvexHull(points)
	ch_vertices = ch.vertices
	return [points[index] for index in ch_vertices]

dimension = 4

all_vertices = pb.get_5_cell_vertices()
# all_vertices = pb.get_orthoplex_vertices(dimension)
# all_vertices = pb.get_cube_vertices(dimension)
# all_vertices = pb.get_24_cell_vertices()
# all_vertices = pb.get_600_cell_vertices()
# all_vertices = pb.get_120_cell_vertices()

random_direction = np.array([0.123, 0.23, 0.76, 0.771])
if len(all_vertices) > 5:
	all_vertices = [v for v in all_vertices if np.inner(v, random_direction) > 0]

# print all_vertices
points = [ np.array(pb.get_origin(dimension)) ]

vertex_count = -1
for v in all_vertices:
	vertex_count += 1
	negative_v = [ - elem for elem in v]
	new_points1 = [[sum(x) for x in zip(p, v)] for p in points]
	new_points2 = [[sum(x) for x in zip(p, negative_v)] for p in points]

	points = np.concatenate((new_points1, new_points2), axis = 0)
	if len(points) > 6 * len(all_vertices):
		points = convex_hull(points)
	# print vertex_count, len(points)

points = convex_hull(points)

scaling = sp.sqrt(1) * 2
points = [p/scaling for p in points]
points.sort(key=lambda x:x[-1] * 1000000 + x[-2] * 10000 + x[-3] * 100 + x[-4])

print "vertex count:", len(points)
print "vertices:"
for p in points:
	print p

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
# 16-cell -> 16 vertices, hypercube
# 8-cell -> 104 vertices, 40 cells: 24 dihedral rhombic dodecahedron + 16 cubes. Can be seen as the rectified 24-cell with eight cubes augmented by pyramids. A rhombic dodecahedron has eight 60-120 rhombi and four squares
# 24-cell -> 192 vertices: truncated 24-cell
# 600-cell -> 7616 (have to use math.sqrt, otherwise too slow)
# 120-cell didn't try
