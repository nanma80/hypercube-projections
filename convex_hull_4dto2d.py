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

high_dimension = 4
low_dimension = 2
phi = (1 + math.sqrt(5.0))/2

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
  print "The number of boundary points of this convex hull: ", len(hull.vertices)
  plt.plot(points[:,0], points[:,1], 'o')
  lines = [[tuple(points[j]) for j in i] for i in edges]
  lc = pltcol.LineCollection(lines)
  plt.axes().add_collection(lc)
  plt.axes().set_aspect('equal')
  plt.show()

def shadow_volume(bases):
  return convex_hull(bases).volume

def negative_volume(bases):
  bases = bases.reshape((low_dimension, high_dimension))
  return -shadow_volume(bases)

def maximize_shadow():
  random_bases = random.rand(low_dimension, high_dimension)
  res_bh = basinhopping(negative_volume, random_bases, disp = True)
  optimal_bases = res_bh.x.reshape((low_dimension, high_dimension))
  orth_optimal_bases = orth(optimal_bases.T).T
  return (- res_bh.fun, orth_optimal_bases)

def get_an_bases(dimension_subspace):
  dimension = dimension_subspace + 1
  base1 = []
  base2 = []
  for index in xrange(dimension):
    base1.append(math.sin(index * 2 * math.pi / dimension))
    base2.append(math.cos(index * 2 * math.pi / dimension))
  return [base1, base2]

# vertices = array(get_5_cell_vertices())
# vertices = array(get_cube_vertices(4))
# vertices = array(get_orthoplex_vertices(4))
vertices = array(get_24_cell_vertices())
# vertices = array(get_120_cell_vertices())
# vertices = array(get_600_cell_vertices())
edges = get_edges(vertices)
print "vertex count:", len(vertices), "edge count:", len(edges)

def main():
  # 6 -> 3
  # known_bases = array([[1, phi, 0, -1, phi, 0], [phi, 0, 1, phi, 0, -1], [0, 1, phi, 0, -1, phi]])
  # print "Volume of the known bases: ", shadow_volume(known_bases)

  # B4
  # t = 1 + math.sqrt(2.0)
  # known_bases = array([[t, t, 1, 1], [1, -1, t, -t]])

  # F4
  a = -1 + math.sqrt(3.0)
  known_bases = array([[1, 1, a, 0], [1, -1, 0, a]])

  # known_bases = array([[1, 0, 0, 0], [0, 1, 0, 0]])

  # H3 approximation
  # known_bases = array([[-0.56386048, -0.43150524, -0.51330528, -0.48206045], [ 0.63692086, -0.73308252,  0.1136261,  -0.20978786]])

  # base1 = [1, 0, 0, 0]
  # base2 = [0, -phi, 1, 0]
  # base3 = [0, 1, phi, 0]
  # known_bases = array([base1, base2])

  # A2 bases for 24-cell
  # a2_bases = get_an_bases(2)
  # known_bases = array([base + [0.0] for base in a2_bases])

  # H4
  # base1 = [(1+math.sqrt(5))*math.sin(math.pi/30), 1, 0, 0]
  # base2 = [0, 0, 2*math.sin(2*math.pi/15), (1+math.sqrt(5))*math.sin(math.pi/15)]
  # base3 = [1, -(1+math.sqrt(5))*math.sin(math.pi/30), 0, 0]
  # base4 = [0, 0, -(1+math.sqrt(5))*math.sin(math.pi/15), 2*math.sin(2*math.pi/15)]
  # known_bases = array([base1, base2])
  # known_bases = array([base3, base4])

  print 'known_bases:'
  print known_bases
  print "Volume of the known bases: ", shadow_volume(known_bases)
  print_convex_hull(known_bases)
  return

  max_shadow, orth_optimal_bases = maximize_shadow()
  print "Volume of max shadow: ", max_shadow
  print "Max achieving bases:"
  print repr(orth_optimal_bases)
  print_convex_hull(orth_optimal_bases)
  inner_products_of_vectors(orth_optimal_bases.T)

if __name__ == '__main__':
  main()

# 120-cell 2D
# area of max shadow:  24.7279388544
# Max achieving bases:
# array([[-0.31146922, -0.25394198,  0.26778851,  0.8756653 ],
#        [-0.07523121, -0.19541302,  0.90863443, -0.36129973]])
# achieved by the first 2 dimensions of H4, or last two
# min shadow 22.2703284592 [10]
# 600-cell 2D
# min shadow 2.7837909157 [10]
# 24-cell min 3.46410208576
# 24-cell max 5.19615242271, achieved by a2 bases
