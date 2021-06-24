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
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

high_dimension = 24
low_dimension = 3
phi = (1 + math.sqrt(5.0))/2
vertex_limit = -1
vertex_limit = 1000

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
  boundary_points = array([points[vector_index]/2.40600382 for vector_index in hull.vertices ])
  print boundary_points
  scaling = 2
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(scaling * points[:,0], scaling * points[:,1], scaling * points[:,2], '.')
  ax.set_aspect('equal')
  ax.axis('off')
  ax.view_init(elev=34, azim=10)
  proj3d.persp_transformation = orthogonal_proj
  plt.show()


def shadow_volume(bases):
  return convex_hull(bases).volume

def negative_volume(bases):
  bases = bases.reshape((low_dimension, high_dimension))
  return - shadow_volume(bases)

def maximize_shadow(known_bases):
  res_bh = basinhopping(negative_volume, known_bases, T=100, niter=1, disp = True)
  optimal_bases = res_bh.x.reshape((low_dimension, high_dimension))
  orth_optimal_bases = orth(optimal_bases.T).T
  return (- res_bh.fun, orth_optimal_bases)

def pad(vectors, target_length):
  for vector_index in range(len(vectors)):
    vectors[vector_index] = vectors[vector_index] + [0] * (target_length - len(vectors[vector_index]))
  return vectors

def get_phi_bases():
  phi = (1 + math.sqrt(5))/2
  base1 = [0, 0, 0, 0, 0, 0, 0, phi, phi, 0, 0, 1, 0, -1]
  base2 = [0, 0, 0, 0, 0, 0, 0, 1, -1, phi, phi, 0, 0, 0]
  base3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, phi, 0, phi]
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
print "All vertex count:", len(all_vertices)
vertices = []
if vertex_limit > 0:
  vertices = sample(all_vertices, vertex_limit)
else:
  vertices = all_vertices
print "Sampled vertex count:", len(vertices)

def main():
  known_bases = array(get_phi_bases())
  # known_bases = array([base + [0.0] for base in a2_bases])
  # known_bases = array(get_bn_bases(4))
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

# phi base
# 393.333194509 all points
# Volume of the known bases:  393.333194509
# The number of boundary points of this convex hull:  42

