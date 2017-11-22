#!/usr/bin/env python

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

# high_dimension = 6
# low_dimension = 3

phi = (1 + math.sqrt(5.0))/2

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

# different initial orientation
def get_24_cell_vertices_2():
  vertices = []
  vertices.extend(get_cube_vertices(4))
  vertices.extend(get_orthoplex_vertices(4))
  for i in xrange(len(vertices)):
    vertices[i] = [math.sqrt(2) * x for x in vertices[i]] # to match lengths of the first version
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

def get_rectified_24_cell_vertices():
  vertices = []
  vertices.extend(get_even_perm_3([1, 1, 2], 0))
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

def projection(vertices, bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  orth_bases = orth(bases.T)
  v2d = np.dot(vertices, orth_bases)
  return v2d

def negative_volume(bases):
  bases = bases.reshape((low_dimension, high_dimension))
  return - shadow_volume(bases)

def maximize_shadow():
  random_bases = random.rand(low_dimension, high_dimension)
  res_bh = basinhopping(negative_volume, random_bases, disp = False)
  optimal_bases = res_bh.x.reshape((low_dimension, high_dimension))
  orth_optimal_bases = orth(optimal_bases.T).T
  return (- res_bh.fun, orth_optimal_bases)


def main():
  # 6 -> 3
  if high_dimension == 6 and low_dimension == 3:
    phi = (1 + math.sqrt(5.0))/2
    known_bases = array([[1, phi, 0, -1, phi, 0], [phi, 0, 1, phi, 0, -1], [0, 1, phi, 0, -1, phi]])
  # print "Volume of the known bases: ", shadow_volume(known_bases)

  # 4 -> 2
  # t = 1 + math.sqrt(2.0)
  # known_bases = array([[t, t, 1, 1], [1, -1, t, -t]])
  # print "Volume of the known bases: ", shadow_volume(known_bases)

  # 8 -> 4
  if high_dimension == 8 and low_dimension == 4:
    known_bases = array([[ 0.33782819, -0.3045361,   0.37740119,  0.48751729,  0.10664892,  0.51613686,   0.2219143,   0.29327104],
                         [ 0.5708417,   0.60604744,  0.15616942, -0.09854069, -0.21650482,  0.22772449,   0.03706754, -0.41550008],
                         [ 0.15099414,  0.11750736, -0.3665696,  -0.46027767,  0.17204068,  0.13761818,   0.64332522,  0.39339532],
                         [-0.19290317,  0.16173212, -0.44587536,  0.20189433, -0.64199191,  0.40349154,  -0.18845184,  0.29426483]])
  vertices = get_cube_vertices(high_dimension)

  projections = np.dot(vertices, known_bases.T)
  ch = ConvexHull(projections)
  projections = [projections[index] for index in ch.vertices]
  convex_hull_vertices = [vertices[index] for index in ch.vertices]

  
  start_index = 119
  # Good indices:
  #   8 8
  # 14 4
  # 26 4
  # 28 4
  # 32 4
  # 39 8
  # 48 4
  # 56 4
  # 71 4
  # 79 4
  # 88 8
  # 95 4
  # 99 4
  # 101 4
  # 113 4
  # 119 8

  v0 = projections[start_index]

  enhanced_projections = []
  for i, v1 in enumerate(projections):
    cos = np.inner(v0, v1)/math.sqrt(np.inner(v0, v0) * np.inner(v1, v1))
    # cos = np.inner(v0, v1)
    enhanced_projections.append((i, cos, v1, convex_hull_vertices[i]))

  enhanced_projections.sort(key = lambda x: x[1])
  threshold = enhanced_projections[-2][1]
  top_vertices = [x for x in enhanced_projections if x[1] > threshold - 0.001 and x[1] < threshold + 0.001]

  for enhanced_projection in enhanced_projections:
    print enhanced_projection

  print len(enhanced_projections)

  print "Volume of the convex hull: ", ch.volume, "Volume of the known bases: ", shadow_volume(known_bases)

  print "For comparison"
  # Seems to be related to the edges of 24-cell (vertices of rectified 24-cell)
  # comparison_vertices = get_rectified_24_cell_vertices()
  # comparison_vertices = get_24_cell_vertices()
  comparison_vertices = [v[2] for v in top_vertices]
  v0 = comparison_vertices[0]

  comparison_enhanced_projections = []
  for i, v1 in enumerate(comparison_vertices):
    cos = np.inner(v0, v1)/math.sqrt(np.inner(v0, v0) * np.inner(v1, v1))
    # cos = np.inner(v0, v1)
    comparison_enhanced_projections.append((i, cos, v1))
  comparison_enhanced_projections.sort(key = lambda x: x[1])

  for x in comparison_enhanced_projections:
    print x


  # max_shadow, orth_optimal_bases = maximize_shadow()
  # print "Volume of max shadow: ", max_shadow
  # print "Max achieving bases:"
  # print orth_optimal_bases
  # inner_products_of_vectors(orth_optimal_bases.T)

if __name__ == '__main__':
  main()
