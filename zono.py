import numpy as np
from scipy.spatial import ConvexHull
import polytope_builder as pb

dimension = 4

# all_vertices = pb.get_5_cell_vertices()
all_vertices = pb.get_orthoplex_vertices(dimension)
# all_vertices = pb.get_cube_vertices(dimension)
# all_vertices = pb.get_24_cell_vertices()
# all_vertices = pb.get_600_cell_vertices()


# todo: randomize the order in all_vertices
# dedup

points = [ np.array(pb.get_origin(dimension)) ]

for v in all_vertices:
	points = points[:]
	new_points = [[sum(x) for x in zip(p, v)] for p in points]
	points = np.concatenate((points, new_points), axis = 0)
	if len(points) > 100:
		ch = ConvexHull(points)
		ch_vertices = ch.vertices
		points = [points[index] for index in ch.vertices]

print len(points)
# for p in points:
# 	print p
