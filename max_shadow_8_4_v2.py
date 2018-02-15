#!/usr/bin/env python
# Special logic for 8D -> 4D to find the best bases in human understandable way

import sys
import itertools
import math
import numpy as np
from scipy import random, array, dot, zeros
from scipy.linalg import orth, det
from scipy.optimize import minimize, basinhopping

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
  a = math.cos(x)/2
  b = math.sin(x)/2
  c = math.cos(y)/2
  d = math.sin(y)/2

  return array([
    [a,  a,  a,  a, b, b, b,  b],
    [b,  b, -b, -b, a, a, -a, -a],
    [c, -c, d, -d, d, -d, c, -c],
    [d, -d, c, -c, -c, c, -d, d]])


def negative_volume(generators):
  return - shadow_volume(get_bases(generators))

def maximize_shadow():
  # generators = [0.378556, 0.326643, 0.100944, 0.489704 ]
  # generators = [0.711917, 0.203285]
  generators = random.rand(2,1)
  # generators = [math.pi/4, math.pi/4]
  print "Initial generator is:"
  print generators
  print "Initial bases are:" 
  print get_bases(generators) * math.sqrt(2) * 2
  print "Initial Volume is:", shadow_volume(get_bases(generators))

  bounds = Bounds([math.pi/4, math.pi/4], [0,0])
  res_bh = basinhopping(negative_volume, generators, disp = False, niter=50, accept_test=bounds)
  print "Optimal generator is:"
  print res_bh.x
  optimal_bases = get_bases(res_bh.x)
  return (- res_bh.fun, optimal_bases)

def main():
  max_shadow, orth_optimal_bases = maximize_shadow()
  print "Volume of max shadow: ", max_shadow

if __name__ == '__main__':
  main()
