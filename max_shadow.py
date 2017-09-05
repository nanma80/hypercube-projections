#!/usr/bin/env python

import sys
import itertools
import math
import numpy as np
from scipy import random, array, dot, zeros
from scipy.linalg import orth, det
from scipy.optimize import minimize, basinhopping

high_dimension = 8
low_dimension = 4

def inner_products_of_vectors(vectors):
  for indices in itertools.combinations(range(len(vectors)), 2):
    print indices, np.inner(vectors[indices[0]], vectors[indices[1]])

def shadow_volume(bases):
  _high_dimension = len(bases.T)
  _low_dimension = len(bases)
  orth_bases = orth(bases.T)
  sum = 0.0
  for indices in itertools.combinations(range(_high_dimension), _low_dimension):
    sub_matrix = orth_bases[indices, :]
    sum += abs(det(sub_matrix))
  return sum

def negative_volume(bases):
  bases = bases.reshape((low_dimension, high_dimension))
  return - shadow_volume(bases)

def maximize_shadow():
  random_bases = random.rand(low_dimension, high_dimension)
  res_bh = basinhopping(negative_volume, random_bases, disp = False)
  optimal_bases = res_bh.x.reshape((low_dimension, high_dimension))
  orth_optimal_bases = orth(optimal_bases.T).T
  return (- res_bh.fun, orth_optimal_bases)


def main():
  # 6 -> 3
  # phi = (1 + math.sqrt(5.0))/2
  # known_bases = array([[1, phi, 0, -1, phi, 0], [phi, 0, 1, phi, 0, -1], [0, 1, phi, 0, -1, phi]])
  # print "Volume of the known bases: ", shadow_volume(known_bases)

  # 4 -> 2
  # t = 1 + math.sqrt(2.0)
  # known_bases = array([[t, t, 1, 1], [1, -1, t, -t]])
  # print "Volume of the known bases: ", shadow_volume(known_bases)
  max_shadow, orth_optimal_bases = maximize_shadow()
  print "Volume of max shadow: ", max_shadow
  print "Max achieving bases:"
  print orth_optimal_bases
  print "Inner products between bases:"
  inner_products_of_vectors(orth_optimal_bases.T)

if __name__ == '__main__':
  main()
