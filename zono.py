import numpy as np
from scipy.spatial import ConvexHull
import polytope_builder as pb

dimension = 4

# all_vertices = pb.get_600_cell_vertices()
all_vertices = pb.get_orthoplex_vertices(dimension)

points = [ pb.get_origin(dimension) ]

for v in all_vertices:
	points = points[:]
	new_points = [[sum(x) for x in zip(p, v)] for p in points]
	points = np.concatenate((points, new_points), axis = 0)
	if len(points) > 100:
		points = ConvexHull(points).points

print len(points)
for p in points:
	print p
