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
# import matplotlib.pyplot as plt
# import matplotlib.collections as pltcol
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d import proj3d


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
    # co = [2 * j - 1 for j in [int(ch) for ch in str][1:]]
    co = [j - 0.5 for j in [int(ch) for ch in str][1:]]
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
  ring2 = get_demicube_vertices(8, True)
  vertices.extend(ring2)
  return vertices

def convex_hull(bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  norm = np.inner(bases[0], bases[0])
  normalizing_factor = math.sqrt(norm)
  # phi = (1 + sqrt(5))/2
  # normalizing_factor = phi
  normalized_bases = bases.T / normalizing_factor

  v2d = np.dot(vertices, normalized_bases)
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

def get_phi_bases():
  phi = (1 + sqrt(5))/2
  bases = [
      [1, phi, 0, -1, phi, 0, 0, 0], 
      [phi, 0, 1, phi, 0, -1, 0, 0], 
      [0, 1, phi, 0, -1, phi, 0, 0],
      [0, 0, 0, 0, 0, 0, phi+1, phi-1]
    ]

  return pad(bases, high_dimension)

def get_optimal_bases(generators):
  x,y = generators
  a = math.cos(x)/2
  b = math.sin(x)/2
  c = math.cos(y)/2
  d = math.sin(y)/2

  return array([
    [a,  a,  a,  a, b, b, b,  b],
    [b,  b, -b, -b, a, a, -a, -a],
    [c, -c, d, -d, d, -d, c, -c],
    [d, -d, c, -c, -c, c, -d, d]])

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

# for dim in [3, 4, 5, 6, 7, 8]:
#   vertices = array(get_demicube_vertices(dim))
#   print dim, len(vertices)
vertices = array(get_cube_vertices(8))
# vertices = array(get_demicube_vertices(8))
# vertices = array(get_4_21_vertices())

# for v in vertices:
#   print repr(v)
print "vertex count: ", len(vertices)
# edges = get_edges(vertices)
# print "edge count: ", len(edges)

# bases = get_e6_bases()
phi_bases = array(get_phi_bases())
# bases = get_bn_bases(5)
# bases = get_an_bases(7)

# print repr(array(phi_bases))


generators = [0.71191743, 0.20328551]
optimal_bases = get_optimal_bases(generators)

# phi = (1 + sqrt(5))/2
# known_bases = known_bases / phi

print "Volume of the known bases: ", shadow_volume(phi_bases)
print "Volume of the optimal bases: ", shadow_volume(optimal_bases)

known_bases = optimal_bases

hull = convex_hull(known_bases)
hull_vertices = [hull.points[vertex_index] for vertex_index in hull.vertices]
hull_vertices.sort(key=lambda x:x[3])

print "# vertices in the convex hull:", len(hull.vertices)
hull_edges = get_edges(hull_vertices)
print "# edges in the convex hull:", len(hull_edges)
for index, vertex in enumerate(hull_vertices):
  print ''.join(["%02d" % (index+1), ': ['] + ["{:10.4f}".format(x) for x in vertex] + [']'])
  # pass

# phi bases:
# vertices in the convex hull: 64
# edges in the convex hull: 120

# optimal bases:
# vertices in ch: 128
# edges of max length: 48. But there must be other edges