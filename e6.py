#!/usr/bin/env python

import datetime
import numpy as np
# import sympy as sp
# from scipy.spatial import ConvexHull
# import polytope_builder as pb
# from sympy import sqrt
from math import sqrt

import sys
import itertools
import math
from scipy import random, array, dot, zeros
from scipy.linalg import orth, det
from scipy.optimize import minimize, basinhopping
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.collections as pltcol
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d


high_dimension = 7
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

def get_edges(vertices):
  target_inner_product = max([np.inner(vertices[0], vertices[i]) for i in xrange(1, len(vertices))])
  print target_inner_product/np.inner(vertices[0], vertices[0])
  edges = []
  for i in xrange(len(vertices)):
    for j in xrange(i+1, len(vertices)):
      inner_prod = np.inner(vertices[i], vertices[j])
      if abs(inner_prod - target_inner_product) < 0.01:
        edges.append([i, j])
  return edges


def get_origin(dimension):
  v = []
  for i in range(dimension):
    v.append(0)
  return v

def get_cube_vertices(dimension):
  vertices = []
  limit = 2 ** dimension
  for i in range(limit):
    str = "{0:b}".format(i + limit)
    co = [2 * j - 1 for j in [int(ch) for ch in str][1:]]
    vertices.append(co)
  return vertices

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

def get_double_non_zero_vertices(dimension):
  vertices = []
  for i in range(dimension):
    for j in range(i+1, dimension):
      for k in range(2):
        for l in range(2):
          vertex = get_origin(dimension)
          vertex[i] = k * 2 - 1
          vertex[j] = l * 2 - 1
          vertices.append(vertex)
  return vertices

def get_demicube_vertices(dimension, alt_mode = False):
  remainder = ((0 if alt_mode else 2) + dimension) % 4
  vertices = [vector for vector in get_cube_vertices(dimension) if (sum(vector) + 8) % 4 == remainder]
  return vertices

def get_2_21_vertices():
  vertices = []
  if high_dimension == 6:
    vertices.append([0, 0, 0, 0, 0, 4 / sqrt(3)])
    ring2 = [vector + [1/sqrt(3)] for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 3]
    vertices.extend(ring2)
    ring3 = [[2 * el for el in vector] + [ - 2 / sqrt(3)] for vector in get_orthoplex_vertices(5)]
    vertices.extend(ring3)
  elif high_dimension == 7:
    ring1 = [[- 6 * el + 1 for el in vector] + [3 * sqrt(2)] for vector in get_orthoplex_vertices(6) if sum(vector) > 0]
    ring2 = [[- 6 * el + 1 for el in vector] + [-3 * sqrt(2)] for vector in get_orthoplex_vertices(6) if sum(vector) > 0]
    vertices.extend(ring1)
    vertices.extend(ring2)
    ring3 = [[ 6 * el - 2 for el in vector] + [0] for vector in get_double_non_zero_vertices(6) if sum(vector) > 1.5]
    vertices.extend(ring3)
  elif high_dimension == 8:
    vertices.append([0, 0, 0, 0, 0] + [4/3.] * 3)
    ring2 = [vector + [1/3.]*3 for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 1]
    vertices.extend(ring2)
    ring3 = [[2 * el for el in vector] + [ - 2/3.]*3 for vector in get_orthoplex_vertices(5)]
    vertices.extend(ring3)
  return vertices

def get_1_22_vertices():
  vertices = []
  if high_dimension == 6:
    ring1 = [[2 * el for el in vector] + [0] for vector in get_double_non_zero_vertices(5)]
    vertices.extend(ring1)
    ring2 = [vector + [sqrt(3)] for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 1]
    vertices.extend(ring2)
    ring3 = [vector + [-sqrt(3)] for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 3]
    vertices.extend(ring3)
  elif high_dimension == 8:
    ring1 = [[2 * el for el in vector] + [0]*3 for vector in get_double_non_zero_vertices(5)]
    vertices.extend(ring1)
    ring2 = [vector + [1]*3 for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 3]
    vertices.extend(ring2)
    ring3 = [vector + [-1]*3 for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 1]
    vertices.extend(ring3)

  return vertices

def get_omega_powers():
  return [
    [1, 0],
    [- 1.0/2.0, sqrt(3)/2.0],
    [- 1.0/2.0, -sqrt(3)/2.0]
  ]

def get_hessian_vertices():
  # The real representation of the Hessian polytope is 2_21
  vertices = []
  for index1 in xrange(3):
    for index2 in xrange(3):
      part1 = [0, 0]
      part2 = get_omega_powers()[index1]
      part3 = [-x for x in get_omega_powers()[index2]]
      vertices.append(part1 + part2 + part3)
      vertices.append(part2 + part3 + part1)
      vertices.append(part3 + part1 + part2)
  return vertices

def convex_hull(bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  orth_bases = orth(bases.T)
  v2d = np.dot(vertices, orth_bases)
  hull = ConvexHull(v2d)
  return hull

def print_convex_hull_2d(bases):
  hull = convex_hull(bases)
  points = hull.points
  plt.plot(points[:,0], points[:,1], 'o')
  lines = [[tuple(points[j]) for j in i] for i in edges]
  lc = pltcol.LineCollection(lines)
  plt.axes().add_collection(lc)
  plt.axes().set_aspect('equal')
  plt.show()

def print_convex_hull_3d(bases):
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
  ax.view_init(elev=90, azim=0)
  proj3d.persp_transformation = orthogonal_proj
  plt.show()

def print_convex_hull(bases):
  if low_dimension == 2:
    print_convex_hull_2d(bases)
  elif low_dimension == 3:
    print_convex_hull_3d(bases)
  else:
    print "low_dimension = ", low_dimension, "not supported"

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

def pad(vectors, target_length):
  for vector_index in range(len(vectors)):
    vectors[vector_index] = vectors[vector_index] + [0] * (target_length - len(vectors[vector_index]))
  return vectors

def get_bases():
  phi = (1 + sqrt(5))/2
  base1 = [phi, phi, 0, 0, -1, 1]
  base2 = [1, -1, phi, -phi, 0, 0]
  base3 = [0, 0, 1, 1, phi, phi]
  return pad([base1, base2, base3], high_dimension)

def get_bases_rearranged():
  phi = (1 + sqrt(5))/2
  base1 = [phi, 0, 1, phi, 0, -1]
  base2 = [0, 1, phi, 0, -1, phi]
  base3 = [1, phi, 0, -1, phi, 0]
  return pad([base1, base2, base3], high_dimension)

def get_e6_bases():
  a = sqrt(3) - 1
  base1 = [1, 1, a, 0, 0, 0]
  base2 = [1, -1, 0, a, 0, 0]
  base3 = [0, 0, 0, 0, 1, 0]
  return pad([base1, base2, base3], high_dimension)

def get_bn_bases(dimension):
  base1 = []
  base2 = []
  for index in xrange(dimension):
    base1.append(math.cos(index * math.pi / dimension + math.pi / dimension / 2))
    base2.append(math.sin(index * math.pi / dimension + math.pi / dimension / 2))
  return pad([base1, base2], high_dimension)

def get_an_bases(dimension_subspace):
  dimension = dimension_subspace + 1
  base1 = []
  base2 = []
  for index in xrange(dimension):
    base1.append(math.sin(index * 2 * math.pi / dimension))
    base2.append(math.cos(index * 2 * math.pi / dimension))
  return pad([base1, base2], high_dimension)

def get_a5_special_bases(): # special A5 for 2_21 and 1_22
  return array([
    [2, 2, 0, 0, 1, sqrt(3)],
    [1,-1, 1, -1, 0, 0]])
  # return array([
  #   [sqrt(3)/2,-sqrt(3)/2, sqrt(3)/2, -sqrt(3)/2, 0, 0],
  #   [1, 1, 0, 0, 1./2, sqrt(3)/2]]) # normalized

# vertices = array(get_cube_vertices(6))
# vertices = array(get_demicube_vertices(6, False))
# vertices = array(get_demicube_vertices(6, True))
# vertices = array(get_1_22_vertices()) # doesn't work well with A5, B6 projections. vertices of the last dimension is probably not standard
vertices = array(get_2_21_vertices()) # doesn't work well with A5, B6 projections. vertices of the last dimension is probably not standard
# vertices = array(get_hessian_vertices()) # equivalent to 2_21?

print "vertex count: ", len(vertices)
# for v in vertices:
#   print repr(v)
edges = get_edges(vertices)
print "edge count: ", len(edges)

# bases = get_e6_bases()
# bases = get_bases()
bases = get_bases_rearranged()
# bases = get_bn_bases(6)
# bases = get_an_bases(5)
# bases = array([base + [0, 0, 0] for base in get_an_bases(2)])
# bases = get_a5_special_bases()

print repr(array(bases))
known_bases = array(bases)
if low_dimension == 2:
  known_bases = array([bases[0], bases[1]])

print "Volume of the known bases: ", shadow_volume(known_bases)
print_convex_hull(known_bases)

# max_shadow, orth_optimal_bases = maximize_shadow()
# print "Volume of max shadow: ", max_shadow
# print "Max achieving bases:"
# print repr(orth_optimal_bases)
# print_convex_hull(orth_optimal_bases)

# 6-cube:
# 64 vertices, 192 edges
# max 2D projection: B6, 14.928
# max 3D projection: get_bases, phi, 1, 0, phi, -1, 0 .... convex hull: rhombic triacontahedron
# 
# 6-demicube:
# 32 vertices, 240 edges
# max 2D projection: A5 original formation (not alt). 13.8564
# max 3D projection: get_bases, phi, 1, 0, phi, -1, 0 .... project from the alt formation. convex hull: icosahedron 30.38666
# 
# 1_22 polytope:
# 72 vertices, 720 edges
# max 2D projection: D4/A2, 20.784
# max 3D projection: D4 projection + another dimension
# convex hull similar to a hexagon prism
# Volume of max shadow:  56.442694636
# 
# 2_21 polytope:
# 27 vertices, 216 edges
# max 2D projection: A5, 10.3923048454
# Max achieving bases: not standard A5, probably because the vertex coordinates are not standard
# [[1, 2, 0, 0,  2, sqrt(3)],
#  [0, 1, 1, 1, -1,       0]]
# max 3D projection: D5 projection + another dimension, 
# convex hull: square pyramid + square antiprism + square pyramid (Gyroelongated square bipyramid)
# Volume of max shadow:  20.0947570825
# 
