#!/usr/bin/env python

import sys
import itertools
import math
import numpy as np
from scipy import random, array, dot, zeros
from scipy.linalg import orth, det
from scipy.optimize import minimize, basinhopping
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.collections as pltcol
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

high_dimension = 4
low_dimension = 3
 
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    # -0.0001 added for numerical stability as suggested in:
    # http://stackoverflow.com/questions/23840756
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,-0.0001,zback]])

def rotate(l, n):
    return l[-n:] + l[:-n]

def get_cube_vertices(dimension):
  vertices = []
  limit = 2 ** dimension
  for i in range(limit):
    str = "{0:b}".format(i + limit)
    co = [1.0 * j - 0.5 for j in [int(ch) for ch in str][1:]]
    vertices.append(co)
  return vertices

def get_even_perm_3(prototype, last):
  prototypes_level1 = []
  limit = 2 ** 3
  for i in range(limit):
    str = "{0:b}".format(i + limit)
    coefficients = [2.0 * j - 1 for j in [int(ch) for ch in str][1:]]
    prototypes_level1.append( [prototype[i] * coefficients[i] for i in xrange(3)] )

  prototypes_level2 = []
  for r in xrange(3):
    for p in prototypes_level1:
      prototypes_level2.append(rotate(p, r) + [last])

  prototypes_level3 = []
  for p in prototypes_level2:
    prototypes_level3.append([p[0], p[1], p[2], p[3]])
    prototypes_level3.append([p[3], p[2], p[1], p[0]])
    prototypes_level3.append([p[2], p[3], p[0], p[1]])
    prototypes_level3.append([p[1], p[0], p[3], p[2]])
  return prototypes_level3

def get_origin(dimension):
  v = []
  for i in range(dimension):
    v.append(0.0)
  return v

def get_orthoplex_vertices(dimension):
  vertices = []
  for i in range(dimension):
    vertex0 = get_origin(dimension)
    vertex0[i] = 1
    vertices.append(vertex0)
    vertex1 = get_origin(dimension)
    vertex1[i] = -1
    vertices.append(vertex1)
  return vertices

def get_24_cell_vertices():
  vertices = []
  for i in range(4):
    for j in range(i+1, 4):
      for k in range(2):
        for l in range(2):
          vertex = get_origin(4)
          vertex[i] = k * 2.0 - 1.0
          vertex[j] = l * 2.0 - 1.0
          vertices.append(vertex)
  return vertices

def get_5_cell_vertices():
  vertices = [
    [1, 1, 1, -1 / math.sqrt(5)],
    [1, -1, -1, -1 / math.sqrt(5)],
    [-1, 1, -1, -1 / math.sqrt(5)],
    [-1, -1, 1, -1 / math.sqrt(5)],
    [0, 0, 0, math.sqrt(5) - 1 / math.sqrt(5)],
  ]
  return vertices

def get_120_cell_vertices():
  vertices = []
  for v in get_24_cell_vertices():
    vertices.append([2 * i for i in v])

  phi = (math.sqrt(5.0) + 1.0) / 2.0
  prototypes = [[1, 1, 1, math.sqrt(5)], [phi, phi, phi, 1/phi/phi], [1/phi, 1/phi, 1/phi, phi*phi]]
  cube_vertices = get_cube_vertices(4)

  for prototype in prototypes:
    prototypes_level1 = []
    for cv in cube_vertices:
      prototypes_level1.append( [prototype[i] * cv[i] * 2 for i in xrange(4)] )
    for p in prototypes_level1:
      for r in xrange(4):
        vertices.append(rotate(p, r))

  vertices.extend(get_even_perm_3([1/phi/phi, 1, phi*phi], 0))
  vertices.extend(get_even_perm_3([1/phi, phi, math.sqrt(5)], 0))
  vertices.extend(get_even_perm_3([1, phi, 2], 1/phi))
  vertices.extend(get_even_perm_3([1, phi, 2], -1/phi))
  return vertices


def get_600_cell_vertices():
  vertices = []
  phi = (math.sqrt(5.0) + 1.0) / 2.0
  prototype = [0.5 * i for i in [phi, 1, 1.0 / phi]]
  vertices.extend(get_even_perm_3(prototype, 0))
  vertices.extend(get_cube_vertices(4))
  vertices.extend(get_orthoplex_vertices(4))
  return vertices

def get_edges(vertices):
  target_inner_product = max([np.inner(vertices[0], vertices[i]) for i in xrange(1, len(vertices))])
  edges = []
  for i in xrange(len(vertices)):
    for j in xrange(i+1, len(vertices)):
      inner_prod = np.inner(vertices[i], vertices[j])
      if abs(inner_prod - target_inner_product) < 0.01:
        edges.append([i, j])
  return edges

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

def convex_hull(bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  orth_bases = orth(bases.T)
  v2d = np.dot(vertices, orth_bases)
  hull = ConvexHull(v2d)
  return hull

def print_convex_hull(bases):
  hull = convex_hull(bases)
  points = hull.points
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(points[:,0], points[:,1], points[:,2], marker = 'o')
  for edge in edges:
    line = [list(points[j]) for j in edge]
    ax.plot([line[0][0], line[1][0]], [line[0][1], line[1][1]], [line[0][2], line[1][2]], color='k')
  ax.set_aspect('equal')
  ax.axis('off')
  proj3d.persp_transformation = orthogonal_proj
  plt.show()

def shadow_volume(bases):
  return convex_hull(bases).volume

def negative_volume(bases):
  bases = bases.reshape((low_dimension, high_dimension))
  return - shadow_volume(bases)

def maximize_shadow():
  random_bases = random.rand(low_dimension, high_dimension)
  res_bh = basinhopping(negative_volume, random_bases, disp = False)
  optimal_bases = res_bh.x.reshape((low_dimension, high_dimension))
  orth_optimal_bases = orth(optimal_bases.T).T
  return (- res_bh.fun, orth_optimal_bases)

def orth_base(bases):
  known_bases = np.append(bases, [[1] * len(bases[0])], axis = 0)
  last_base = np.linalg.solve(known_bases, array([0, 0, 0, 1]))
  # last_base = last_base / math.sqrt(np.inner( last_base, last_base)) * math.sqrt(np.inner(known_bases[0], known_bases[0]))
  last_base = last_base / max([abs(i) for i in last_base])
  return last_base

# vertices = array(get_5_cell_vertices()) # edge first
# volume: 3.26598632371
# [1, -1, -1, 3/sqrt(5)]

# vertices = array(get_cube_vertices(4)) # vertex first
# volume: 2.0  [1, 1, 1, 1]

# vertices = array(get_orthoplex_vertices(4)) # vertex first
# 1.3333 [1, 0, 0, 0]

vertices = array(get_24_cell_vertices()) # what?
# 24-cell: target volume: 7.05533682951 vector [-0.1889823  -0.18898229  0.18898219  0.94491117]
# 0.94491117/0.18898219 = 5. So the vector is [-1, -1, 1, 5]

# vertices = array(get_120_cell_vertices()) # unclear
# Volume of max shadow:  87.3688309937
# [ 0.14818048 -0.23976104 -0.13253656  1.        ]
# [ 1.          0.79944109  0.05013976  0.05013967]
# [ 0.18033988  1.         -0.78521826  0.61803399]
# [-0.04095612  1.          0.04095609 -0.30602924]

# vertices = array(get_600_cell_vertices()) # close to vertex first (3.55 vs 3.53)
# Volume of max shadow:  3.55713925244
# [ 0.30444186  1.         -0.57012138  0.12543673]
# [ -1, 0, 0.795320722, 0.491535219] and note that 0.795320722 = phi * 0.491535219

edges = get_edges(vertices)
# print "vertex count:", len(vertices), "edge count:", len(edges)

def main():
  # trivial bases
  # known_bases = array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
  # known_bases = array([[1, 1, 1, 1], [1, -1, -1, 1], [-1, 1, -1, 1]])

  # B4
  # t = 1 + math.sqrt(2.0)
  # known_bases = array([[t, t, 1, 1], [-1, 1, t, -t], [1, 1, -t, -t]])

  # F4
  # a = -1 + math.sqrt(3.0)
  # known_bases = array([[1, 1, a, 0], [1, -1, 0, a], [a, 0, -1, -1]])

  # best base for 24-cell
  b = 1./3 # or 2, 5, 1/3, generating the same volume
  known_bases = array([[b, 1, 1, 1], [-1, b, 1, -1], [-1, -1, b, 1]])

  print "Volume of the known bases: ", shadow_volume(known_bases)
  # print_convex_hull(known_bases)
  # return

  max_shadow, orth_optimal_bases = maximize_shadow()
  print "Volume of max shadow: ", max_shadow
  # print "Max achieving bases:"
  # print orth_optimal_bases
  print "orth vector to the optimal bases:"
  print orth_base(orth_optimal_bases)
  print_convex_hull(orth_optimal_bases)
  # inner_products_of_vectors(orth_optimal_bases.T)


if __name__ == '__main__':
  # for i in xrange(10):
  main()
