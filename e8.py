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


high_dimension = 8
low_dimension = 4

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

def get_4_21_vertices():
  vertices = []
  ring1 = [[2 * el for el in vector] for vector in get_double_non_zero_vertices(8)]
  vertices.extend(ring1)
  ring2 = get_demicube_vertices(8)
  vertices.extend(ring2)
  return vertices

def convex_hull(bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  orth_bases = orth(bases.T)
  v2d = np.dot(vertices, orth_bases)
  hull = ConvexHull(v2d)
  return hull

def dist(x, y):
  return np.sqrt(np.sum((x-y)**2))

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
  bases = [
      [1, phi, 0, -1, phi, 0, 0, 0], 
      [phi, 0, 1, phi, 0, -1, 0, 0], 
      [0, 1, phi, 0, -1, phi, 0, 0],
      [0, 0, 0, 0, 0, 0, phi+1, phi-1]
    ]

  return pad(bases, high_dimension)

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

for dim in [3, 4, 5, 6, 7, 8]:
  vertices = array(get_demicube_vertices(dim))
  print dim, len(vertices)
# vertices = array(get_cube_vertices(8))
# vertices = array(get_demicube_vertices(8))
vertices = array(get_4_21_vertices())

# for v in vertices:
#   print repr(v)
print "vertex count: ", len(vertices)
edges = get_edges(vertices)
print "edge count: ", len(edges)

# bases = get_e6_bases()
bases = get_bases()
# bases = get_bn_bases(5)
# bases = get_an_bases(7)

print repr(array(bases))
known_bases = array(bases)
if low_dimension == 2:
  known_bases = array([bases[0], bases[1]])
if low_dimension == 3:
  known_bases = array([bases[0], bases[1], bases[2]])

print "Volume of the known bases: ", shadow_volume(known_bases)
print_convex_hull(known_bases)
ch = convex_hull(known_bases)
points_on_convex_hull = [ch.points[index] for index in ch.vertices]
edge_length = min([dist(points_on_convex_hull[0], points_on_convex_hull[i]) for i in xrange(1, len(points_on_convex_hull))])
# for v in points_on_convex_hull:
#   print repr(v)

print edge_length
print 26.475425 * (edge_length ** 4), ch.volume

# max_shadow, orth_optimal_bases = maximize_shadow()
# print "Volume of max shadow: ", max_shadow
# print "Max achieving bases:"
# print repr(orth_optimal_bases)
# print_convex_hull(orth_optimal_bases)


# 4_21:
# bases based on phi: Volume of the known bases:  129.4427191
# optimal bases: Volume of max shadow:  142.08103671
# Max achieving bases:
# array([[-0.20428321, -0.59858509,  0.06389811, -0.59506347,  0.24978869,
#          0.05710992,  0.24540529, -0.34044246],
#        [ 0.42419816, -0.46386136, -0.06555465, -0.18770753, -0.15416367,
#         -0.46957864, -0.41506279,  0.38575818],
#        [ 0.55479698, -0.11375784, -0.11333372, -0.08582672, -0.06081787,
#          0.60916705,  0.40573804,  0.3458931 ],
#        [ 0.27089508,  0.12824633,  0.76264977,  0.02096803,  0.55432917,
#         -0.08949402,  0.00405637,  0.11308381]])
