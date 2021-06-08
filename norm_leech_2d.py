#!/usr/bin/env python

import sys
import itertools
import math
import csv
import numpy as np
from random import sample
from scipy import random, array, dot, zeros
from scipy.linalg import orth, det
from scipy.optimize import minimize, basinhopping
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.collections as pltcol

high_dimension = 24
low_dimension = 2
phi = (1 + math.sqrt(5.0))/2
vertex_limit = -1
vertex_limit = 10000

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

def convex_hull(bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  orth_bases = orth(bases.T)
  v2d = np.dot(vertices, orth_bases)
  hull = ConvexHull(v2d)
  return hull

def plot_points(bases):
  hull = convex_hull(bases)
  points = hull.points
  print "The number of boundary points of this convex hull: ", len(hull.vertices)
  plt.plot(points[:,0], points[:,1], '.', markersize=2)
  plt.axes().set_aspect('equal')
  plt.show()

def max_squared_norm(bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  orth_bases = orth(bases.T)
  v2d = np.dot(vertices, orth_bases)
  return max(np.dot(v, v) for v in v2d)

def objective_function(bases):
  bases = bases.reshape((low_dimension, high_dimension))
  return max_squared_norm(bases)

def minimize_max_squared_norm(known_bases):
  # random_bases = random.rand(low_dimension, high_dimension)
  res_bh = basinhopping(objective_function, known_bases, T=100, niter=10, disp = True)
  optimal_bases = res_bh.x.reshape((low_dimension, high_dimension))
  orth_optimal_bases = orth(optimal_bases.T).T
  return (res_bh.fun, orth_optimal_bases)

def pad(vectors, target_length):
  for vector_index in range(len(vectors)):
    vectors[vector_index] = vectors[vector_index] + [0] * (target_length - len(vectors[vector_index]))
  return vectors

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

def get_leech_vertices():
  data_file_name = 'external_data/minvects.csv'
  vertices = []
  line_number = 0
  with open(data_file_name) as csvfile:
    reader = csv.reader(csvfile, delimiter = ' ')
    for row in reader:
      row_numbers = [int(i) for i in row]
      # if max([abs(i) for i in row_numbers]) == 4:
      vertices.append(row_numbers)
      line_number += 1
      # if vertex_limit > 0 and line_number >= vertex_limit:
        # break
  return vertices

all_vertices = get_leech_vertices()
vertices = []
if vertex_limit > 0:
  vertices = sample(all_vertices, vertex_limit)
else:
  vertices = all_vertices
print "vertex count:", len(vertices)

def main():
  # known_bases = array([base + [0.0] for base in a2_bases])
  known_bases = array(get_bn_bases(4))
  known_bases = array([[-0.10653864,  0.14554017, -0.24932734,  0.28836727,  0.1792456 ,
         0.22652015, -0.22292096, -0.16912753, -0.06721476,  0.25347954,
        -0.05447918,  0.22474247, -0.27849967,  0.11673808,  0.28642218,
         0.09259233,  0.08631201, -0.26843499, -0.12039588, -0.01877529,
         0.29097719,  0.28572895,  0.24813508,  0.18473565],
       [ 0.27081771, -0.25296869, -0.13883825, -0.00046725,  0.22453279,
         0.17577163, -0.18709428, -0.21857623, -0.29127737, -0.13411774,
         0.28956325, -0.19204695,  0.05854167, -0.2680833 , -0.03329871,
        -0.27292654,  0.26502643, -0.1212536 ,  0.26485656, -0.2856895 ,
        -0.00900591, -0.02694499, -0.1288563 ,  0.22522294]])
  print(known_bases)

  print "max squared norm of the 2d points from the known bases: ", max_squared_norm(known_bases)
  # plot_points(known_bases)

  min_max_sq_norm, orth_optimal_bases = minimize_max_squared_norm(known_bases)
  print "minimized max squared norm of the 2d points: ", min_max_sq_norm
  print "Min achieving bases:"
  print repr(orth_optimal_bases)
  plot_points(orth_optimal_bases)

if __name__ == '__main__':
  main()

# random 1000: 7.13
# 1104 vertices of the third class: 
# minimized max squared norm of the 2d points:  5.44456094283
# Min achieving bases:
# array([[-0.10653864,  0.14554017, -0.24932734,  0.28836727,  0.1792456 ,
#          0.22652015, -0.22292096, -0.16912753, -0.06721476,  0.25347954,
#         -0.05447918,  0.22474247, -0.27849967,  0.11673808,  0.28642218,
#          0.09259233,  0.08631201, -0.26843499, -0.12039588, -0.01877529,
#          0.29097719,  0.28572895,  0.24813508,  0.18473565],
#        [ 0.27081771, -0.25296869, -0.13883825, -0.00046725,  0.22453279,
#          0.17577163, -0.18709428, -0.21857623, -0.29127737, -0.13411774,
#          0.28956325, -0.19204695,  0.05854167, -0.2680833 , -0.03329871,
#         -0.27292654,  0.26502643, -0.1212536 ,  0.26485656, -0.2856895 ,
#         -0.00900591, -0.02694499, -0.1288563 ,  0.22522294]])
# The number of boundary points of this convex hull:  44
