import numpy as np
from scipy.spatial import ConvexHull
import polytope_builder as pb

dimension = 4

# all_vertices = pb.get_5_cell_vertices()
# all_vertices = pb.get_orthoplex_vertices(dimension)
all_vertices = pb.get_cube_vertices(dimension)
# all_vertices = pb.get_24_cell_vertices()
# all_vertices = pb.get_600_cell_vertices()
# all_vertices = pb.get_120_cell_vertices()

random_direction = np.array([0.123, 0.23, 0.76, 0.771])
all_vertices = [v for v in all_vertices if np.inner(v, random_direction) > 0]

points = [ np.array(pb.get_origin(dimension)) ]

vertex_count = -1
for v in all_vertices:
	negative_v = [ - elem for elem in v]
	points = points[:]
	new_points1 = [[sum(x) for x in zip(p, v)] for p in points]
	new_points2 = [[sum(x) for x in zip(p, negative_v)] for p in points]

	points = np.concatenate((new_points1, new_points2), axis = 0)
	if len(points) > 100:
		ch = ConvexHull(points)
		ch_vertices = ch.vertices
		points = [points[index] for index in ch.vertices]
	print vertex_count, len(points)
	vertex_count += 1


print len(points)
# for p in points:
# 	print p


# 5-cell -> 8
# 16-cell -> 16 (hypercube)
# 8-cell -> 104?
# 24-cell -> 192
# 600-cell -> 7616 (have to use math.sqrt, otherwise too slow)
# 120-cell didn't try
