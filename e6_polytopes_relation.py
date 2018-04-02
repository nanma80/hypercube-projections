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


high_dimension = 6
low_dimension = 2

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

def get_6_demicube_vertices():
  vertices = [vector for vector in get_cube_vertices(6) if (sum(vector) + 8) % 4 == 0]
  return vertices

def get_6_demicube_vertices_alt():
  vertices = [vector for vector in get_cube_vertices(6) if (sum(vector) + 8) % 4 == 2]
  return vertices

def get_2_21_vertices():
  vertices = []
  if high_dimension == 6:
    vertices.append([0, 0, 0, 0, 0, 4 / sqrt(3)])
    ring2 = [vector + [1/sqrt(3)] for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 3]
    vertices.extend(ring2)
    ring3 = [[2 * el for el in vector] + [ - 2 / sqrt(3)] for vector in get_orthoplex_vertices(5)]
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

def dist(x,y):   
  return np.sqrt(np.sum((x-y)**2))

# vertices = array(get_cube_vertices(6))
# vertices = array(get_6_demicube_vertices())
# vertices = array(get_6_demicube_vertices_alt())
vertices122 = array(get_1_22_vertices()) # doesn't work well with A5, B6 projections. vertices of the last dimension is probably not standard
vertices221 = [[3./4 * el for el in vertex] for vertex in array(get_2_21_vertices())] # doesn't work well with A5, B6 projections. vertices of the last dimension is probably not standard
vertices221_opposite = [[-3./4 * el for el in vertex] for vertex in array(get_2_21_vertices())]
vertices221 = vertices221 + vertices221_opposite

print "vertex count: ", len(vertices122)
print "vertex count: ", len(vertices221)
edges122 = get_edges(vertices122)
edges221 = get_edges(vertices221)
print "edge count: ", len(edges122)
print "edge count: ", len(edges221)

target_dist = min([dist(vertices221[0], v) for v in vertices122])
group_a, group_b = vertices221, vertices122
# group_b, group_a = vertices221, vertices122

for va in group_a:
  count = 0
  for vb in group_b:
    distance = dist(va, vb)
    if abs(distance - target_dist) < 0.01:
      count += 1
  print count, va

