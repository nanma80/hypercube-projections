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
vertex_limit = 1000

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

def print_convex_hull(bases):
  hull = convex_hull(bases)
  points = hull.points
  print "The number of boundary points of this convex hull: ", len(hull.vertices)
  plt.plot(points[:,0], points[:,1], '.', markersize=2)
  plt.axes().set_aspect('equal')
  plt.show()

def shadow_volume(bases):
  return convex_hull(bases).volume

def negative_volume(bases):
  bases = bases.reshape((low_dimension, high_dimension))
  return - shadow_volume(bases)

def maximize_shadow(known_bases):
  res_bh = basinhopping(negative_volume, known_bases, T=100, niter=100, disp = True)
  optimal_bases = res_bh.x.reshape((low_dimension, high_dimension))
  orth_optimal_bases = orth(optimal_bases.T).T
  return (- res_bh.fun, orth_optimal_bases)

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
  # known_bases = array([[-1.51203363e-07, -5.77757193e-08, -1.60925410e-08,
  #        4.11199837e-01,  1.38878440e-07, -6.57779631e-09,
  #       -1.28925099e-07, -7.63799721e-08, -5.75251995e-01,
  #        1.24866730e-07,  4.58985868e-08, -6.97526581e-01,
  #       -2.40936051e-07,  1.61515362e-07,  8.87248787e-08,
  #       -5.15954223e-08, -1.16002176e-01,  7.35137736e-08,
  #       -4.26356098e-08,  2.47091607e-08,  1.92259735e-08,
  #        3.18078251e-07,  2.27489040e-07,  1.22825777e-07],
  #      [-1.27949344e-07, -1.51720427e-07, -1.91433519e-08,
  #       -5.75251905e-01,  2.41147732e-08, -7.67538455e-08,
  #        1.81344750e-08,  2.33090848e-08, -4.11199979e-01,
  #       -1.08195805e-07, -1.51765524e-07,  1.16002340e-01,
  #       -3.89924410e-08,  1.57861467e-08, -4.13338160e-09,
  #        1.21350020e-07, -6.97526544e-01, -5.41514110e-08,
  #       -9.84793153e-08, -7.24993736e-08,  3.36016628e-08,
  #       -8.64660982e-08,  4.44426533e-08, -3.73727783e-08]])

  print(known_bases)

  print "Volume of the known bases: ", shadow_volume(known_bases)
  print_convex_hull(known_bases)
  # return

  max_shadow, orth_optimal_bases = maximize_shadow(known_bases)
  print "Volume of max shadow: ", max_shadow
  print "Max achieving bases:"
  print repr(orth_optimal_bases)
  print_convex_hull(orth_optimal_bases)

if __name__ == '__main__':
  main()

# first 1000 points: convex hull 8-gon. 72.86122023
# array([[-1.47067514e-01,  5.79418777e-08, -5.96694704e-08,
#         -6.66574998e-09,  5.23930002e-08, -7.52689266e-08,
#         -1.47067733e-01,  7.57595376e-08, -2.51064674e-01,
#         -2.51064895e-01, -4.98008741e-01, -4.26537573e-01,
#          6.36796947e-08,  1.09586344e-07,  9.28602910e-08,
#         -2.79469981e-01,  7.03006068e-08, -3.50941339e-01,
#          6.99062478e-08, -6.76194312e-08, -2.79470207e-01,
#         -2.79469982e-01, -1.47067698e-01, -1.47067675e-01],
#        [ 1.38888237e-01, -1.14302650e-07, -5.96801806e-08,
#         -2.11086208e-08,  1.11822060e-07, -6.00882031e-08,
#          1.38888436e-01, -9.56444258e-08,  4.64202655e-01,
#          4.64202464e-01,  1.67511433e-01, -2.30693577e-01,
#          1.92434212e-07, -4.34514972e-08,  1.08098309e-07,
#         -3.69582208e-01,  6.56441966e-08,  2.86230911e-02,
#          7.61119894e-08, -8.84807057e-08, -3.69582277e-01,
#         -3.69582251e-01,  1.38888460e-01,  1.38888542e-01]])
# first 10000 points: convex hull 8-gon. 77.2548339959
# last 1104 points: max(abs) == 4: 77.2548
# array([[-1.51203363e-07, -5.77757193e-08, -1.60925410e-08,
#          4.11199837e-01,  1.38878440e-07, -6.57779631e-09,
#         -1.28925099e-07, -7.63799721e-08, -5.75251995e-01,
#          1.24866730e-07,  4.58985868e-08, -6.97526581e-01,
#         -2.40936051e-07,  1.61515362e-07,  8.87248787e-08,
#         -5.15954223e-08, -1.16002176e-01,  7.35137736e-08,
#         -4.26356098e-08,  2.47091607e-08,  1.92259735e-08,
#          3.18078251e-07,  2.27489040e-07,  1.22825777e-07],
#        [-1.27949344e-07, -1.51720427e-07, -1.91433519e-08,
#         -5.75251905e-01,  2.41147732e-08, -7.67538455e-08,
#          1.81344750e-08,  2.33090848e-08, -4.11199979e-01,
#         -1.08195805e-07, -1.51765524e-07,  1.16002340e-01,
#         -3.89924410e-08,  1.57861467e-08, -4.13338160e-09,
#          1.21350020e-07, -6.97526544e-01, -5.41514110e-08,
#         -9.84793153e-08, -7.24993736e-08,  3.36016628e-08,
#         -8.64660982e-08,  4.44426533e-08, -3.73727783e-08]])
# Only 4 vectors, angle between them: pi/4
# nothing but array(get_bn_bases(4))
# array([[0,0],
#        [0,0],
#        [0,0],
#        [0.290762193, -0.406764523],
#        [0,0],
#        [0,0],
#        [0,0],
#        [0,0],
#        [-0.406764587, -0.290762294],
#        [0,0],
#        [0,0],
#        [-0.493225775,  0.0820260412],
#        [0,0],
#        [0,0],
#        [0,0],
#        [0,0],
#        [-0.0820259253, -0.493225749],
#        [0,0],
#        [0,0],
#        [0,0],
#        [0,0],
#        [0,0],
#        [0,0],
#        [0,0]])
