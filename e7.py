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

def get_2_31_vertices():
  vertices = []
  if high_dimension == 7:
    vertices.append([0, 0, 0, 0, 0, 0, 2 * sqrt(2)])
    vertices.append([0, 0, 0, 0, 0, 0, - 2 * sqrt(2)])
    ring1 = [[2 * el for el in vector] + [0] for vector in get_double_non_zero_vertices(6)]
    vertices.extend(ring1)
    ring2 = [vector + [sqrt(2)] for vector in get_demicube_vertices(6, True)]
    vertices.extend(ring2)
    ring3 = [vector + [-sqrt(2)] for vector in get_demicube_vertices(6, True)]
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

# vertices = array(get_cube_vertices(7))
vertices = array(get_demicube_vertices(7, False))
# vertices = array(get_demicube_vertices(7, True))
# vertices = array(get_2_31_vertices())

print "vertex count: ", len(vertices)
# for v in vertices:
#   print repr(v)
edges = get_edges(vertices)
print "edge count: ", len(edges)

# bases = get_e6_bases()
bases = get_bases()
# bases = get_bases_rearranged()
# bases = get_bn_bases(7)
# bases = get_an_bases(6)
# bases = array([base + [0, 0, 0] for base in get_an_bases(2)])
# bases = get_a5_special_bases()

print repr(array(bases))
known_bases = array(bases)
if low_dimension == 2:
  known_bases = array([bases[0], bases[1]])

print "Volume of the known bases: ", shadow_volume(known_bases)
print_convex_hull(known_bases)

max_shadow, orth_optimal_bases = maximize_shadow()
print "Volume of max shadow: ", max_shadow
print "Max achieving bases:"
print repr(orth_optimal_bases)
print_convex_hull(orth_optimal_bases)

# 7-cube:
# Projection to 2D: B7/A6, Volume of max shadow:  17.5251450701
# Projection to 3D: 
# Volume of max shadow:  44.8310806031
# Max achieving bases:
# array([[-0.42186143, -0.31731367,  0.11278828,  0.00672709,  0.49525002,
#          0.64509627,  0.21715603],
#        [-0.36990147,  0.26681348, -0.11609854, -0.63633067,  0.27102376,
#         -0.11121148, -0.53643826],
#        [ 0.34138303, -0.51379829, -0.63650359,  0.10145105,  0.34111507,
#         -0.08572404, -0.28342848]])

# 7-demicube:
# Projection to 2D: B7/A6, Volume of max shadow:  15.7896101137
# Projection to 3D: 
# Volume of max shadow:  39.1981082676
# Max achieving bases:
# array([[-0.53319699, -0.14035373, -0.02392288,  0.45808228,  0.55648033,
#         -0.41798626,  0.03474544],
#        [-0.2396595 , -0.60797219, -0.48020182, -0.40420859,  0.19432319,
#          0.36157606,  0.10225468],
#        [ 0.28599597, -0.19284761,  0.44916599,  0.21548184,  0.2980004 ,
#          0.3605986 ,  0.64342677]])
# 
# 2_31 polytope:
# 126 vertices, 2016 edges
# When projected by get_bases, the convex hull is an icosidodecahedron. 
# It contains less points than the projection of 4_21 
# (not icosidodecahedron, but the canonical projection of 600-cell)
# Projection to 2D:
# Volume of max shadow:  20.7846096908
# D4/B3/A2/G2
# Projection to 3D:
# Volume of max shadow:  56.442694636
# Max achieving bases:
# array([[-0.43995093, -0.38346406, -0.35522075, -0.16798342, -0.38346409,
#         -0.25271383, -0.54230001],
#        [-0.52059141,  0.11384073,  0.43105663,  0.33704509,  0.11384064,
#         -0.61460286,  0.16099472],
#        [ 0.24942623,  0.09817574,  0.02255056, -0.77411475,  0.09817611,
#         -0.54723891,  0.13884137]])
