#!/usr/bin/env python
# Special logic for 8D -> 4D to find the best bases in human understandable way
# compare to v2, we changed the definition of the projection matrix, so that the two solutions are symmetric

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

class Bounds(object):
  def __init__(self, xmax=[1,1], xmin=[0,0] ):
    self.xmax = np.array(xmax)
    self.xmin = np.array(xmin)
  def __call__(self, **kwargs):
    x = kwargs["x_new"]
    tmax = bool(np.all(x <= self.xmax))
    tmin = bool(np.all(x >= self.xmin))
    return tmax and tmin

def add_columns(bases):
  new_bases = np.zeros((low_dimension, high_dimension))
  new_bases[0, high_dimension-2] = math.sqrt(0.5)
  new_bases[1, high_dimension-1] = math.sqrt(0.5)
  new_bases[:, :-2] = bases
  return new_bases

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

def get_bases(generators):
  x,y = generators
  cx = math.cos(x)
  sx = math.sin(x)
  cy = math.cos(y)
  sy = math.sin(y)

  return array([[sx, sx, cx, cx, -sx, -sx, cx, cx],
    [cx, cx, sx, sx, cx, cx, -sx, -sx],
    [cy, -cy, sy, -sy, sy, -sy, cy, -cy],
    [sy, -sy, cy, -cy, -cy, cy, -sy, sy]])

def negative_volume(generators):
  return - shadow_volume(get_bases(generators))

def maximize_shadow():
  # generators = [0.378556, 0.326643, 0.100944, 0.489704 ]
  # generators = [0.711917, 0.203285]
  generators = random.rand(2,1)
  # generators = [math.pi/4, math.pi/4]
  # print "Initial generator is:"
  # print generators
  # print "Initial bases are:" 
  # print get_bases(generators) * math.sqrt(2) * 2
  # print "Initial Volume is:", shadow_volume(get_bases(generators))

  bounds = Bounds([0.08, math.pi/8], [0,0])
  # bounds = Bounds([180, 180], [0,0])
  res_bh = basinhopping(negative_volume, generators, disp = False, niter=100, accept_test=bounds)
  a1 = res_bh.x[0]
  a2 = res_bh.x[1]
  if a1 > 0.08 or a2 > math.pi/8:
    raise ValueError("basinhopping gets a solution out of the bounds")
  print "Optimal generator is:"
  print [a1, a2]
  print "in degrees, a1, a2, sum, diff"
  print [180/math.pi * a for a in [a1, a2, a1+a2, a1-a2]]

  optimal_bases = get_bases(res_bh.x)
  return (- res_bh.fun, optimal_bases)

def convex_hull(vertices, bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  norm = np.inner(bases[0], bases[0])
  normalizing_factor = math.sqrt(norm)
  # phi = (1 + sqrt(5))/2
  # normalizing_factor = phi
  normalized_bases = bases.T #/ normalizing_factor

  v2d = np.dot(vertices, normalized_bases)
  hull = ConvexHull(v2d)
  return hull

def get_cube_vertices(dimension):
  vertices = []
  limit = 2 ** dimension
  for i in range(limit):
    str = "{0:b}".format(i + limit)
    # co = [2 * j - 1 for j in [int(ch) for ch in str][1:]]
    co = [j - 0.5 for j in [int(ch) for ch in str][1:]]
    vertices.append(co)
  return vertices

def main():
  max_shadow, orth_optimal_bases = maximize_shadow()
  print "Volume of max shadow: ", max_shadow
  print "target shadow volume: ", 7.84468782041
  # print "optimal base"
  # print orth_optimal_bases

  vertices = array(get_cube_vertices(8))
  hull = convex_hull(vertices, orth_optimal_bases)
  hull_vertices = [hull.points[vertex_index] for vertex_index in hull.vertices]
  hull_vertices.sort(key=lambda x:x[3])

  print "# vertices in the convex hull:", len(hull.vertices)
  for index, vertex in enumerate(hull_vertices):
    print ''.join(["%02d" % (index+1), ': ['] + ["{:10.4f}".format(x) for x in vertex] + [']'])
  # print hull.vertices
  
if __name__ == '__main__':
  main()

# in degrees, a1, a2, sum, diff
# [11.647401655103769, 4.210135109188659, 15.857536764292428, 7.437266545915109]
# [4.210136353618086, 11.64740236003748, 15.857538713655567, -7.437266006419395]
