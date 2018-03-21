#!/usr/bin/env python

import datetime
import numpy as np
# import sympy as sp
# from scipy.spatial import ConvexHull
# import polytope_builder as pb
# from sympy import sqrt
from math import sqrt

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


def get_2_21_vertices():
  vertices = []
  vertices.append([0, 0, 0, 0, 0, 4 / sqrt(3)])
  ring2 = [vector + [1/sqrt(3)] for vector in get_cube_vertices(5) if sum(vector) % 4 == 1]
  vertices.extend(ring2)

  ring3 = [[2 * el for el in vector] + [ - 2 / sqrt(3)] for vector in get_orthoplex_vertices(5)]
  vertices.extend(ring3)
  return vertices

def get_1_22_vertices():
  vertices = []
  ring1 = [[2 * el for el in vector] + [0] for vector in get_double_non_zero_vertices(5)]
  vertices.extend(ring1)

  ring2 = [vector + [sqrt(3)] for vector in get_cube_vertices(5) if sum(vector) % 4 == 3]
  vertices.extend(ring2)
  ring3 = [vector + [-sqrt(3)] for vector in get_cube_vertices(5) if sum(vector) % 4 == 1]
  vertices.extend(ring3)
  return vertices

# vertices = get_2_21_vertices()
vertices = get_1_22_vertices()

print len(vertices)

edges = get_edges(vertices)
print len(edges)
