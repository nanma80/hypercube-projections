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

def negative_volume(bases):
  bases = bases.reshape((low_dimension, high_dimension))
  return - shadow_volume(bases)

def maximize_shadow():
  random_bases = array([[ 0.33782819, -0.3045361,   0.37740119,  0.48751729,  0.10664892,  0.51613686,
   0.2219143,   0.29327104],
 [ 0.5708417,   0.60604744,  0.15616942, -0.09854069, -0.21650482,  0.22772449,
   0.03706754, -0.41550008],
 [ 0.15099414,  0.11750736, -0.3665696,  -0.46027767,  0.17204068,  0.13761818,
   0.64332522, 0.39339532],
 [-0.19290317,  0.16173212, -0.44587536,  0.20189433, -0.64199191,  0.40349154,
  -0.18845184,  0.29426483]])
  res_bh = minimize(negative_volume, random_bases, method='Nelder-Mead')
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
  inner_products_of_vectors(orth_optimal_bases.T)

if __name__ == '__main__':
  main()
