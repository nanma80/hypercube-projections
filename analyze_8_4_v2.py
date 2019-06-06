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

def inner_products_of_vectors(vectors):
  print "Length of the generators:"
  for v in vectors:
    print math.sqrt(np.inner(v, v))

  cos_tuples = []
  print "Cosines between the generators:"
  for indices in itertools.combinations(range(len(vectors)), 2):
    v0 = vectors[indices[0]]
    v1 = vectors[indices[1]]
    cos = np.inner(v0, v1)/math.sqrt(np.inner(v0, v0) * np.inner(v1, v1))
    cos_tuples.append((indices, cos, math.acos(cos) / math.pi * 180))
  cos_tuples.sort(key=lambda x:x[2])
  for t in cos_tuples:
    print t

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

def main():
  generators = [0.71191743, 0.20328551]
  bases = get_bases(generators)
  x,y = generators
  a = math.cos(x)
  b = math.sin(x)
  c = math.cos(y)
  d = math.sin(y)
  print generators
  print [g / math.pi * 180 for g in generators]

  # max_shadow, orth_bases = maximize_shadow()
  # print "Volume of max shadow: ", max_shadow
  print "orth_bases: "
  print bases

  inner_products_of_vectors(bases.T)


if __name__ == '__main__':
  main()
