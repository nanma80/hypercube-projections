#!/usr/bin/env python
# Special logic for 8D -> 4D to find the best bases in human understandable way
# compare to v2, we changed the definition of the projection matrix, so that the two solutions are symmetric

import sys
import itertools
import math
import numpy as np
from scipy import random, array, dot, zeros
from scipy.linalg import orth, det
from scipy.optimize import minimize, basinhopping
from scipy.spatial import ConvexHull

high_dimension = 8
low_dimension = 4
epsilon = 0.0000001;

class Bounds(object):
  def __init__(self, xmax=[1,1], xmin=[0,0] ):
    self.xmax = np.array(xmax)
    self.xmin = np.array(xmin)
  def __call__(self, **kwargs):
    x = kwargs["x_new"]
    tmax = bool(np.all(x <= self.xmax))
    tmin = bool(np.all(x >= self.xmin))
    return tmax and tmin

def add_columns(bases):
  new_bases = np.zeros((low_dimension, high_dimension))
  new_bases[0, high_dimension-2] = math.sqrt(0.5)
  new_bases[1, high_dimension-1] = math.sqrt(0.5)
  new_bases[:, :-2] = bases
  return new_bases

def inner_products_of_vectors(vectors):
  print "Length of the generators:"
  for v in vectors:
    print math.sqrt(np.inner(v, v))

  print "Cosines between the generators:"
  for indices in itertools.combinations(range(len(vectors)), 2):
    v0 = vectors[indices[0]]
    v1 = vectors[indices[1]]
    cos = np.inner(v0, v1)/math.sqrt(np.inner(v0, v0) * np.inner(v1, v1))
    print indices, cos

def shadow_volume(bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  orth_bases = orth(bases.T)
  sum = 0.0
  for indices in itertools.combinations(range(_high_dimension), _low_dimension):
    sub_matrix = orth_bases[indices, :]
    sum += abs(det(sub_matrix))
  return sum

def get_bases(generators):
  x,y = generators
  cx = math.cos(x)
  sx = math.sin(x)
  cy = math.cos(y)
  sy = math.sin(y)

  return array(
    [
      [sx,  sx, cx,  cx, -sx, -sx,  cx,  cx],
      [cx,  cx, sx,  sx,  cx,  cx, -sx, -sx],
      [sy, -sy, cy, -cy, -cy,  cy,  sy, -sy],
      [cy, -cy, sy, -sy,  sy, -sy, -cy,  cy]
    ])

def negative_volume(generators):
  return - shadow_volume(get_bases(generators))

def maximize_shadow():
  # generators = [0.378556, 0.326643, 0.100944, 0.489704 ]
  # generators = [0.711917, 0.203285]
  generators = random.rand(2,1) * 0.08
  # generators = [math.pi/4, math.pi/4]
  # print "Initial generator is:"
  # print generators
  # print "Initial bases are:" 
  # print get_bases(generators) * math.sqrt(2) * 2
  # print "Initial Volume is:", shadow_volume(get_bases(generators))

  bounds = Bounds([0.08, math.pi/8], [0,0])
  # bounds = Bounds([180, 180], [0,0])
  res_bh = basinhopping(negative_volume, generators, disp = False, niter=100, accept_test=bounds)
  a1 = res_bh.x[0]
  a2 = res_bh.x[1]
  if a1 > 0.08 or a2 > math.pi/8:
    raise ValueError("basinhopping gets a solution out of the bounds")
  print "Optimal generator is:"
  print [a1, a2]
  print "in degrees, a1, a2, sum, diff"
  print [180/math.pi * a for a in [a1, a2, a1+a2, a1-a2]]

  optimal_bases = get_bases(res_bh.x)
  return (- res_bh.fun, optimal_bases)

def convex_hull(vertices, bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  norm = np.inner(bases[0], bases[0])
  normalizing_factor = math.sqrt(norm)
  # phi = (1 + sqrt(5))/2
  # normalizing_factor = phi
  normalized_bases = bases.T / normalizing_factor

  v2d = np.dot(vertices, normalized_bases)
  hull = ConvexHull(v2d)
  return hull, v2d

def get_cube_vertices(dimension):
  vertices = []
  limit = 2 ** dimension
  for i in range(limit):
    str = "{0:b}".format(i + limit)
    # co = [2 * j - 1 for j in [int(ch) for ch in str][1:]]
    co = [j - 0.5 for j in [int(ch) for ch in str][1:]]
    vertices.append(co)
  return vertices

def main():
  max_shadow, orth_optimal_bases = maximize_shadow()
  print "Volume of max shadow: ", max_shadow
  print "target shadow volume: ", 7.84468782041
  # print "optimal base"
  # print orth_optimal_bases

  vertices = array(get_cube_vertices(8))
  hull, projected_vertices = convex_hull(vertices, orth_optimal_bases)
  print "convex hull volume: ", hull.volume
  hull_vertices = [hull.points[vertex_index] for vertex_index in hull.vertices]
  hull_vertices.sort(key=lambda x:x[3])

  print "# vertices in the convex hull:", len(hull.vertices)
  # for index, vertex in enumerate(hull_vertices):
    # print ''.join(["%02d" % (index+1), ': ['] + ["{:10.4f}".format(x) for x in vertex] + [']'])
  # print hull.vertices

  # before dedup, > 800 equations
  # print hull.equations
  # print len(hull.equations)

  rounded_equations = array([[np.round(x, 15) for x in v] for v in hull.equations])
  # rounded_equations = [tuple(row) for row in rounded_equations]
  unique_equations = np.unique(rounded_equations, axis=0)

  # print unique_equations
  print len(unique_equations)

  # offsets = [np.round(v[4], 6) for v in unique_equations]
  # print sorted(np.unique(offsets))

  # after dedup, 112 cells, 128 vertices. Each cell has 8 vertices. Each vertex belong to 7 cells? no. Some vertex belong to more. some fewer.

  ids_in_cells = []

  for equation in unique_equations:
    ids_in_cell = [vid for vid in range(len(projected_vertices)) if dot( np.concatenate((projected_vertices[vid], array([1]))) , equation ) > - epsilon]
    ids_in_cells.append(ids_in_cell)
    # print ids_in_cell

  # print ids_in_cells
  # each cell has 8 vertices. Each cell has central symmetry.
  cell_centers = []
  for ids in ids_in_cells:
    projected_vertices_in_cell = [projected_vertices[id] for id in ids]
    center = array(np.average(projected_vertices_in_cell, axis=0))
    cell_centers.append(center)
    # relative_vertices = array([ np.round(v - center, 10) for v in projected_vertices_in_cell])
    # negative_relative_vertices = array([ np.round(- v + center, 10) for v in projected_vertices_in_cell])
    # all_relative_vertices = np.concatenate((relative_vertices, negative_relative_vertices))
    # print(len(all_relative_vertices))
    # print( len(np.unique(all_relative_vertices, axis=0)) )
    # all 8. all cells are cross polytopes
    # print( np.unique(all_relative_vertices, axis=0) )

    # print(center)
    # print(array([ v - center for v in projected_vertices_in_cell] ))

    # print(array([ np.linalg.norm(v- center) for v in projected_vertices_in_cell] ))

    # break

  cell_id_groups = [[0, 1, 2, 3, 44, 45, 46, 47, 64, 65, 66, 67, 108, 109, 110, 111], [4, 5, 6, 7, 48, 49, 50, 51, 60, 61, 62, 63, 104, 105, 106, 107], [8, 9, 10, 11, 28, 29, 30, 31, 80, 81, 82, 83, 100, 101, 102, 103], [12, 13, 14, 15, 24, 25, 26, 27, 84, 85, 86, 87, 96, 97, 98, 99], [16, 17, 18, 19, 40, 41, 42, 43, 68, 69, 70, 71, 92, 93, 94, 95], [20, 21, 22, 23, 36, 37, 38, 39, 72, 73, 74, 75, 88, 89, 90, 91], [32, 33, 34, 35, 52, 53, 54, 55, 56, 57, 58, 59, 76, 77, 78, 79]]

  for group in cell_id_groups:
    cell_centers_in_group = [ cell_centers[id] for id in group ]
    group_hull = ConvexHull(cell_centers_in_group)
    # rounded_equations = array([[np.round(x, 15) for x in v] for v in group_hull.equations])
    # rounded_equations = [tuple(row) for row in rounded_equations]
    # unique_equations = np.unique(rounded_equations, axis=0)

    # print unique_equations
    # print group_hull.equations
    # print group_hull.points
    # print 'length of equations in the group is', len(group_hull.equations)
    cells = []
    faces = []
    edges = []
    vertices = []
    face_ids = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    edge_ids = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]

    for equation in group_hull.equations:
      # for point in group_hull.points:
        # print(point)
        # gap = array([dot( np.concatenate((point, array([1]))), equation) ])
        # print(gap)
        # return
      cell = [point_id for point_id in range(len(group_hull.points)) if abs(dot( np.concatenate((group_hull.points[point_id], array([1]))), equation)) < epsilon ]
      # print(cell)
      cells.append(cell)
      faces.extend( [[ cell[id] for id in ids ] for ids in face_ids ])
      edges.extend( [[ cell[id] for id in ids ] for ids in edge_ids ])
      vertices.extend(cell)
    
    unique_faces = np.unique(faces, axis = 0)
    unique_edges = np.unique(edges, axis = 0)
    unique_vertices = np.unique(vertices)
    # print sorted([ len([eid for eid in range(len(edges)) if np.array_equal(array(edges[eid]), unique_edge) ]) for unique_edge in unique_edges])
    vertices_in_edges = []
    for edge in unique_edges:
      vertices_in_edges.extend(edge)
    print sorted([ len([vertex for vertex in vertices_in_edges if vertex == unique_vertex ]) for unique_vertex in unique_vertices])
    # print sorted([ len([vertex for vertex in vertices if vertex == unique_vertex ]) for unique_vertex in unique_vertices])

    # print([
    #   (len(cells)),
    #   (len(faces)),
    #   (len(unique_faces)),
    #   (len(edges)),
    #   (len(unique_edges)),
    #   (len(vertices)),
    #   (len(unique_vertices))
    # ])

    # (len(faces))
    # (len(np.unique(faces, axis = 0)))

    # print(len(edges))
    # print(len(np.unique(edges, axis = 0)))

    # return


  # for projected_vertex in projected_vertices:
  #   equations_on_vertex = [eid for eid in range(len(unique_equations)) if dot( np.concatenate( (array(projected_vertex), array([1])) ), unique_equations[eid] ) > - epsilon]
  #   print(len(equations_on_vertex))
    # break
  # print gap.max()


  return

  

  # print hull.simplices
  # print len(hull.simplices)

  gap = array([dot( np.concatenate((projected_vertex, array([1]))) , hull.equations[100]) for projected_vertex in projected_vertices])
  # print gap.max()

  

if __name__ == '__main__':
  main()

# in degrees, a1, a2, sum, diff
# [11.647401655103769, 4.210135109188659, 15.857536764292428, 7.437266545915109]
# [4.210136353618086, 11.64740236003748, 15.857538713655567, -7.437266006419395]

# Optimal generator is:
# [0.07348072421223484, 0.20328550404073642]
# in degrees, a1, a2, sum, diff
# [4.210135372925817, 11.64740141772384, 15.857536790649657, -7.437266044798022]
# Volume of max shadow:  7.84468782041
# target shadow volume:  7.84468782041
# convex hull volume:    7.84468782041

