GlowScript 2.7 VPython

phi = (1 + sqrt(5.0))/2
high_dimension = 8

def initialize():
    scene.fov = 0.001 # simulate orthographic projection

def inner(a, b):
    if len(a) != len(b):
        raise ValueError('Lengths do not match: ' + a + ' ' + b)
    sum = 0.0
    for i in range(len(a)):
        sum += 1.0 * a[i] * b[i]
    return sum

def rotate(l, n):
    return l[-n:] + l[:-n]

def get_cube_vertices(dimension):
    vertices = []
    limit = 2 ** dimension
    for i in range(limit):
        str = "{0:b}".format(i + limit)
        co = [2 * j - 1 for j in [ch for ch in str][1:]]
        vertices.append(co)
    return vertices

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

def project_to_3d(vertices, bases):
    v3d = []
    for v in vertices:
        p = [inner(v, base) for base in bases]
        v3d.append(vec(p[0], p[1], p[2]))
    return v3d

def get_edges(vertices, target_inner_product):
    edges = []
    for i in range(len(vertices)):
        for j in range(i+1, len(vertices)):
            inner_prod = inner(vertices[i], vertices[j])
            if abs(inner_prod - target_inner_product) < 0.01:
                edges.append([i, j])
    return edges

def draw_wireframe(vectors, edges):
    initialize()
    vertex_size = 0.2
    edge_size = vertex_size / 3
    for v in vectors:
        sphere(pos = v, radius = vertex_size)
    
    for edge in edges:
        cylinder(pos = vectors[edge[0]], axis = vectors[edge[1]] - vectors[edge[0]], radius = edge_size)

def draw_wireframe(vectors, edges, vertex_is_visible):
    initialize()
    vertex_size = 0.4
    edge_size = vertex_size / 3
    for v in vectors:
        sphere(pos = v, radius = vertex_size)
    
    for edge in edges:
        cylinder(pos = vectors[edge[0]], axis = vectors[edge[1]] - vectors[edge[0]], radius = edge_size)


def get_edges(vertices):
  target_inner_product = max([inner(vertices[0], vertices[i]) for i in range(1, len(vertices))])
  edges = []
  for i in range(len(vertices)):
    for j in range(i+1, len(vertices)):
      inner_prod = inner(vertices[i], vertices[j])
      if abs(inner_prod - target_inner_product) < 0.01:
        edges.append([i, j])
  return edges

def pad(vectors, target_length):
  for vector_index in range(len(vectors)):
    vectors[vector_index] = vectors[vector_index] + [0] * (target_length - len(vectors[vector_index]))
  return vectors

def get_bases():
  # on an icosahedron, take one vertex and its 5 neighbors
  phi = (1 + sqrt(5))/2
  base1 = [phi, phi, 0,   0,   -1,   1]
  base2 = [1,   -1,  phi, -phi, 0,   0]
  base3 = [0,    0,  1,    1,   phi, phi]
  # base1 = [1,   0,   phi, phi, 0,   -1]
  # base2 = [0,   phi, 1,   -1,  phi, 0]
  # base3 = [phi, 1,   0,    0,  -1,  phi]
  return pad([base1, base2, base3], high_dimension)

def get_bases_v2():
  # project the v2 421 vertices into the most symmetric projection
  phi = (1 + sqrt(5))/2
  # base1 = [1, 2*phi, -1, 1, 0, 1]
  # base2 = [phi-1, 0, phi+1, -phi+1, 0, phi+1]
  # base3 = [phi, 0, phi, phi, 2, -phi]

  a = 2 * phi
  b = 2
  base1 = [0, 0, a, 0, 0, b]
  base2 = [0, a, 0, b, 0, 0]
  base3 = [a, 0, 0, 0, b, 0]

  return pad([base1, base2, base3], high_dimension)


def get_bases_rearranged():
  # common way to project 6-cube
  phi = (1 + sqrt(5))/2
  base1 = [phi, 0, 1, phi, 0, -1]
  base2 = [0, 1, phi, 0, -1, phi]
  base3 = [1, phi, 0, -1, phi, 0]
  return pad([base1, base2, base3], high_dimension)

def get_e6_bases():
   a = sqrt(3) - 1
   base1 = [1, 1, a, 0, 0, 0]
   base2 = [1, -1, 0, a, 0, 0]
   base3 = [0, 0, 0, 0, 0, 0]
   return [base1, base2, base3]

def get_2_21_bases():
  return [
    [-0.268831,   -0.61440537,  0.37838379,  0.5570654,   0.01253407, -0.31077755],
    [-0.1992443,   0.55622989,  0.41359146,  0.34388504,  0.5614607 ,  0.21530695],
    [-0.76856797, -0.12561368,  0.18676016, -0.41295667, -0.18362164,  0.3929313 ]
    ]

def get_1_22_bases():
  return [
  [  4.84989250e-01,  -1.66784395e-08,   8.48707194e-01,   4.73343185e-02, 2.05526129e-01,  -5.66027314e-08],
  [ -6.41530225e-01,   3.05536559e-08,   1.93551007e-01,   4.06746492e-01, 6.20914059e-01,  -4.53142517e-08],
  [ -4.58658012e-01,   9.89070742e-08,   3.15241078e-01,  -8.30337029e-01, -2.82189321e-02,   5.47028903e-08]]

def get_6_demicube_vertices():
  vertices = [vector for vector in get_cube_vertices(6) if (sum(vector) + 8) % 4 == 0]
  return vertices

def get_6_demicube_vertices_alt():
  vertices = [vector for vector in get_cube_vertices(6) if (sum(vector) + 8) % 4 == 2]
  return vertices

def get_demicube_vertices(dimension, alt_mode = False):
  remainder = ((0 if alt_mode else 2) + dimension) % 4
  vertices = [vector for vector in get_cube_vertices(dimension) if (sum(vector) + 8) % 4 == remainder]
  return vertices

def get_2_21_vertices():
  vertices = []
  if high_dimension == 6:
    vertices.append([0, 0, 0, 0, 0, 4 / sqrt(3)])
    ring2 = [vector + [1/sqrt(3)] for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 3]
    vertices.extend(ring2)
    ring3 = [[2 * el for el in vector] + [ - 2 / sqrt(3)] for vector in get_orthoplex_vertices(5)]
    vertices.extend(ring3)
  elif high_dimension == 8:
    vertices.append([0, 0, 0, 0, 0] + [4/3.] * 3)
    ring2 = [vector + [1/3.]*3 for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 1]
    vertices.extend(ring2)
    ring3 = [[2 * el for el in vector] + [ - 2/3.]*3 for vector in get_orthoplex_vertices(5)]
    vertices.extend(ring3)
  return vertices

def get_1_22_vertices():
  vertices = []
  if high_dimension == 6:
    ring1 = [[2 * el for el in vector] + [0] for vector in get_double_non_zero_vertices(5)]
    vertices.extend(ring1)
    ring2 = [vector + [sqrt(3)] for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 1]
    vertices.extend(ring2)
    ring3 = [vector + [-sqrt(3)] for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 3]
    vertices.extend(ring3)
  elif high_dimension == 8:
    ring1 = [[2 * el for el in vector] + [0]*3 for vector in get_double_non_zero_vertices(5)]
    vertices.extend(ring1)
    ring2 = [vector + [1]*3 for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 3]
    vertices.extend(ring2)
    ring3 = [vector + [-1]*3 for vector in get_cube_vertices(5) if (sum(vector) + 8) % 4 == 1]
    vertices.extend(ring3)

  return vertices

def get_4_21_vertices():
  vertices = []
  ring1 = [[2 * el for el in vector] for vector in get_double_non_zero_vertices(8)]
  vertices.extend(ring1)
  ring2 = get_demicube_vertices(8, True)
  vertices.extend(ring2)
  return vertices

def get_4_21_vertices_v2():
  # from integral octonions
  return [
    [-2, 0, 0, 0, 0, 0, 0, 0],
    [-1, -1, -1, 0, 0, 0, 0, -1],
    [-1, -1, -1, 0, 0, 0, 0, 1],
    [-1, -1, 0, -1, 0, 0, -1, 0],
    [-1, -1, 0, -1, 0, 0, 1, 0],
    [-1, -1, 0, 0, -1, -1, 0, 0],
    [-1, -1, 0, 0, -1, 1, 0, 0],
    [-1, -1, 0, 0, 1, -1, 0, 0],
    [-1, -1, 0, 0, 1, 1, 0, 0],
    [-1, -1, 0, 1, 0, 0, -1, 0],
    [-1, -1, 0, 1, 0, 0, 1, 0],
    [-1, -1, 1, 0, 0, 0, 0, -1],
    [-1, -1, 1, 0, 0, 0, 0, 1],
    [-1, 0, -1, -1, 0, -1, 0, 0],
    [-1, 0, -1, -1, 0, 1, 0, 0],
    [-1, 0, -1, 0, -1, 0, -1, 0],
    [-1, 0, -1, 0, -1, 0, 1, 0],
    [-1, 0, -1, 0, 1, 0, -1, 0],
    [-1, 0, -1, 0, 1, 0, 1, 0],
    [-1, 0, -1, 1, 0, -1, 0, 0],
    [-1, 0, -1, 1, 0, 1, 0, 0],
    [-1, 0, 0, -1, -1, 0, 0, -1],
    [-1, 0, 0, -1, -1, 0, 0, 1],
    [-1, 0, 0, -1, 1, 0, 0, -1],
    [-1, 0, 0, -1, 1, 0, 0, 1],
    [-1, 0, 0, 0, 0, -1, -1, -1],
    [-1, 0, 0, 0, 0, -1, -1, 1],
    [-1, 0, 0, 0, 0, -1, 1, -1],
    [-1, 0, 0, 0, 0, -1, 1, 1],
    [-1, 0, 0, 0, 0, 1, -1, -1],
    [-1, 0, 0, 0, 0, 1, -1, 1],
    [-1, 0, 0, 0, 0, 1, 1, -1],
    [-1, 0, 0, 0, 0, 1, 1, 1],
    [-1, 0, 0, 1, -1, 0, 0, -1],
    [-1, 0, 0, 1, -1, 0, 0, 1],
    [-1, 0, 0, 1, 1, 0, 0, -1],
    [-1, 0, 0, 1, 1, 0, 0, 1],
    [-1, 0, 1, -1, 0, -1, 0, 0],
    [-1, 0, 1, -1, 0, 1, 0, 0],
    [-1, 0, 1, 0, -1, 0, -1, 0],
    [-1, 0, 1, 0, -1, 0, 1, 0],
    [-1, 0, 1, 0, 1, 0, -1, 0],
    [-1, 0, 1, 0, 1, 0, 1, 0],
    [-1, 0, 1, 1, 0, -1, 0, 0],
    [-1, 0, 1, 1, 0, 1, 0, 0],
    [-1, 1, -1, 0, 0, 0, 0, -1],
    [-1, 1, -1, 0, 0, 0, 0, 1],
    [-1, 1, 0, -1, 0, 0, -1, 0],
    [-1, 1, 0, -1, 0, 0, 1, 0],
    [-1, 1, 0, 0, -1, -1, 0, 0],
    [-1, 1, 0, 0, -1, 1, 0, 0],
    [-1, 1, 0, 0, 1, -1, 0, 0],
    [-1, 1, 0, 0, 1, 1, 0, 0],
    [-1, 1, 0, 1, 0, 0, -1, 0],
    [-1, 1, 0, 1, 0, 0, 1, 0],
    [-1, 1, 1, 0, 0, 0, 0, -1],
    [-1, 1, 1, 0, 0, 0, 0, 1],
    [0, -2, 0, 0, 0, 0, 0, 0],
    [0, -1, -1, -1, -1, 0, 0, 0],
    [0, -1, -1, -1, 1, 0, 0, 0],
    [0, -1, -1, 0, 0, -1, -1, 0],
    [0, -1, -1, 0, 0, -1, 1, 0],
    [0, -1, -1, 0, 0, 1, -1, 0],
    [0, -1, -1, 0, 0, 1, 1, 0],
    [0, -1, -1, 1, -1, 0, 0, 0],
    [0, -1, -1, 1, 1, 0, 0, 0],
    [0, -1, 0, -1, 0, -1, 0, -1],
    [0, -1, 0, -1, 0, -1, 0, 1],
    [0, -1, 0, -1, 0, 1, 0, -1],
    [0, -1, 0, -1, 0, 1, 0, 1],
    [0, -1, 0, 0, -1, 0, -1, -1],
    [0, -1, 0, 0, -1, 0, -1, 1],
    [0, -1, 0, 0, -1, 0, 1, -1],
    [0, -1, 0, 0, -1, 0, 1, 1],
    [0, -1, 0, 0, 1, 0, -1, -1],
    [0, -1, 0, 0, 1, 0, -1, 1],
    [0, -1, 0, 0, 1, 0, 1, -1],
    [0, -1, 0, 0, 1, 0, 1, 1],
    [0, -1, 0, 1, 0, -1, 0, -1],
    [0, -1, 0, 1, 0, -1, 0, 1],
    [0, -1, 0, 1, 0, 1, 0, -1],
    [0, -1, 0, 1, 0, 1, 0, 1],
    [0, -1, 1, -1, -1, 0, 0, 0],
    [0, -1, 1, -1, 1, 0, 0, 0],
    [0, -1, 1, 0, 0, -1, -1, 0],
    [0, -1, 1, 0, 0, -1, 1, 0],
    [0, -1, 1, 0, 0, 1, -1, 0],
    [0, -1, 1, 0, 0, 1, 1, 0],
    [0, -1, 1, 1, -1, 0, 0, 0],
    [0, -1, 1, 1, 1, 0, 0, 0],
    [0, 0, -2, 0, 0, 0, 0, 0],
    [0, 0, -1, -1, 0, 0, -1, -1],
    [0, 0, -1, -1, 0, 0, -1, 1],
    [0, 0, -1, -1, 0, 0, 1, -1],
    [0, 0, -1, -1, 0, 0, 1, 1],
    [0, 0, -1, 0, -1, -1, 0, -1],
    [0, 0, -1, 0, -1, -1, 0, 1],
    [0, 0, -1, 0, -1, 1, 0, -1],
    [0, 0, -1, 0, -1, 1, 0, 1],
    [0, 0, -1, 0, 1, -1, 0, -1],
    [0, 0, -1, 0, 1, -1, 0, 1],
    [0, 0, -1, 0, 1, 1, 0, -1],
    [0, 0, -1, 0, 1, 1, 0, 1],
    [0, 0, -1, 1, 0, 0, -1, -1],
    [0, 0, -1, 1, 0, 0, -1, 1],
    [0, 0, -1, 1, 0, 0, 1, -1],
    [0, 0, -1, 1, 0, 0, 1, 1],
    [0, 0, 0, -2, 0, 0, 0, 0],
    [0, 0, 0, -1, -1, -1, -1, 0],
    [0, 0, 0, -1, -1, -1, 1, 0],
    [0, 0, 0, -1, -1, 1, -1, 0],
    [0, 0, 0, -1, -1, 1, 1, 0],
    [0, 0, 0, -1, 1, -1, -1, 0],
    [0, 0, 0, -1, 1, -1, 1, 0],
    [0, 0, 0, -1, 1, 1, -1, 0],
    [0, 0, 0, -1, 1, 1, 1, 0],
    [0, 0, 0, 0, -2, 0, 0, 0],
    [0, 0, 0, 0, 0, -2, 0, 0],
    [0, 0, 0, 0, 0, 0, -2, 0],
    [0, 0, 0, 0, 0, 0, 0, -2],
    [0, 0, 0, 0, 0, 0, 0, 2],
    [0, 0, 0, 0, 0, 0, 2, 0],
    [0, 0, 0, 0, 0, 2, 0, 0],
    [0, 0, 0, 0, 2, 0, 0, 0],
    [0, 0, 0, 1, -1, -1, -1, 0],
    [0, 0, 0, 1, -1, -1, 1, 0],
    [0, 0, 0, 1, -1, 1, -1, 0],
    [0, 0, 0, 1, -1, 1, 1, 0],
    [0, 0, 0, 1, 1, -1, -1, 0],
    [0, 0, 0, 1, 1, -1, 1, 0],
    [0, 0, 0, 1, 1, 1, -1, 0],
    [0, 0, 0, 1, 1, 1, 1, 0],
    [0, 0, 0, 2, 0, 0, 0, 0],
    [0, 0, 1, -1, 0, 0, -1, -1],
    [0, 0, 1, -1, 0, 0, -1, 1],
    [0, 0, 1, -1, 0, 0, 1, -1],
    [0, 0, 1, -1, 0, 0, 1, 1],
    [0, 0, 1, 0, -1, -1, 0, -1],
    [0, 0, 1, 0, -1, -1, 0, 1],
    [0, 0, 1, 0, -1, 1, 0, -1],
    [0, 0, 1, 0, -1, 1, 0, 1],
    [0, 0, 1, 0, 1, -1, 0, -1],
    [0, 0, 1, 0, 1, -1, 0, 1],
    [0, 0, 1, 0, 1, 1, 0, -1],
    [0, 0, 1, 0, 1, 1, 0, 1],
    [0, 0, 1, 1, 0, 0, -1, -1],
    [0, 0, 1, 1, 0, 0, -1, 1],
    [0, 0, 1, 1, 0, 0, 1, -1],
    [0, 0, 1, 1, 0, 0, 1, 1],
    [0, 0, 2, 0, 0, 0, 0, 0],
    [0, 1, -1, -1, -1, 0, 0, 0],
    [0, 1, -1, -1, 1, 0, 0, 0],
    [0, 1, -1, 0, 0, -1, -1, 0],
    [0, 1, -1, 0, 0, -1, 1, 0],
    [0, 1, -1, 0, 0, 1, -1, 0],
    [0, 1, -1, 0, 0, 1, 1, 0],
    [0, 1, -1, 1, -1, 0, 0, 0],
    [0, 1, -1, 1, 1, 0, 0, 0],
    [0, 1, 0, -1, 0, -1, 0, -1],
    [0, 1, 0, -1, 0, -1, 0, 1],
    [0, 1, 0, -1, 0, 1, 0, -1],
    [0, 1, 0, -1, 0, 1, 0, 1],
    [0, 1, 0, 0, -1, 0, -1, -1],
    [0, 1, 0, 0, -1, 0, -1, 1],
    [0, 1, 0, 0, -1, 0, 1, -1],
    [0, 1, 0, 0, -1, 0, 1, 1],
    [0, 1, 0, 0, 1, 0, -1, -1],
    [0, 1, 0, 0, 1, 0, -1, 1],
    [0, 1, 0, 0, 1, 0, 1, -1],
    [0, 1, 0, 0, 1, 0, 1, 1],
    [0, 1, 0, 1, 0, -1, 0, -1],
    [0, 1, 0, 1, 0, -1, 0, 1],
    [0, 1, 0, 1, 0, 1, 0, -1],
    [0, 1, 0, 1, 0, 1, 0, 1],
    [0, 1, 1, -1, -1, 0, 0, 0],
    [0, 1, 1, -1, 1, 0, 0, 0],
    [0, 1, 1, 0, 0, -1, -1, 0],
    [0, 1, 1, 0, 0, -1, 1, 0],
    [0, 1, 1, 0, 0, 1, -1, 0],
    [0, 1, 1, 0, 0, 1, 1, 0],
    [0, 1, 1, 1, -1, 0, 0, 0],
    [0, 1, 1, 1, 1, 0, 0, 0],
    [0, 2, 0, 0, 0, 0, 0, 0],
    [1, -1, -1, 0, 0, 0, 0, -1],
    [1, -1, -1, 0, 0, 0, 0, 1],
    [1, -1, 0, -1, 0, 0, -1, 0],
    [1, -1, 0, -1, 0, 0, 1, 0],
    [1, -1, 0, 0, -1, -1, 0, 0],
    [1, -1, 0, 0, -1, 1, 0, 0],
    [1, -1, 0, 0, 1, -1, 0, 0],
    [1, -1, 0, 0, 1, 1, 0, 0],
    [1, -1, 0, 1, 0, 0, -1, 0],
    [1, -1, 0, 1, 0, 0, 1, 0],
    [1, -1, 1, 0, 0, 0, 0, -1],
    [1, -1, 1, 0, 0, 0, 0, 1],
    [1, 0, -1, -1, 0, -1, 0, 0],
    [1, 0, -1, -1, 0, 1, 0, 0],
    [1, 0, -1, 0, -1, 0, -1, 0],
    [1, 0, -1, 0, -1, 0, 1, 0],
    [1, 0, -1, 0, 1, 0, -1, 0],
    [1, 0, -1, 0, 1, 0, 1, 0],
    [1, 0, -1, 1, 0, -1, 0, 0],
    [1, 0, -1, 1, 0, 1, 0, 0],
    [1, 0, 0, -1, -1, 0, 0, -1],
    [1, 0, 0, -1, -1, 0, 0, 1],
    [1, 0, 0, -1, 1, 0, 0, -1],
    [1, 0, 0, -1, 1, 0, 0, 1],
    [1, 0, 0, 0, 0, -1, -1, -1],
    [1, 0, 0, 0, 0, -1, -1, 1],
    [1, 0, 0, 0, 0, -1, 1, -1],
    [1, 0, 0, 0, 0, -1, 1, 1],
    [1, 0, 0, 0, 0, 1, -1, -1],
    [1, 0, 0, 0, 0, 1, -1, 1],
    [1, 0, 0, 0, 0, 1, 1, -1],
    [1, 0, 0, 0, 0, 1, 1, 1],
    [1, 0, 0, 1, -1, 0, 0, -1],
    [1, 0, 0, 1, -1, 0, 0, 1],
    [1, 0, 0, 1, 1, 0, 0, -1],
    [1, 0, 0, 1, 1, 0, 0, 1],
    [1, 0, 1, -1, 0, -1, 0, 0],
    [1, 0, 1, -1, 0, 1, 0, 0],
    [1, 0, 1, 0, -1, 0, -1, 0],
    [1, 0, 1, 0, -1, 0, 1, 0],
    [1, 0, 1, 0, 1, 0, -1, 0],
    [1, 0, 1, 0, 1, 0, 1, 0],
    [1, 0, 1, 1, 0, -1, 0, 0],
    [1, 0, 1, 1, 0, 1, 0, 0],
    [1, 1, -1, 0, 0, 0, 0, -1],
    [1, 1, -1, 0, 0, 0, 0, 1],
    [1, 1, 0, -1, 0, 0, -1, 0],
    [1, 1, 0, -1, 0, 0, 1, 0],
    [1, 1, 0, 0, -1, -1, 0, 0],
    [1, 1, 0, 0, -1, 1, 0, 0],
    [1, 1, 0, 0, 1, -1, 0, 0],
    [1, 1, 0, 0, 1, 1, 0, 0],
    [1, 1, 0, 1, 0, 0, -1, 0],
    [1, 1, 0, 1, 0, 0, 1, 0],
    [1, 1, 1, 0, 0, 0, 0, -1],
    [1, 1, 1, 0, 0, 0, 0, 1],
    [2, 0, 0, 0, 0, 0, 0, 0]
]

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

def get_3_21_vertices(alt_mode = False):
  vertices = []
  if high_dimension == 7:
    ring1 = [[2 * el for el in vector] + [sqrt(2)] for vector in get_orthoplex_vertices(6)]
    vertices.extend(ring1)
    ring2 = [[2 * el for el in vector] + [-sqrt(2)] for vector in get_orthoplex_vertices(6)]
    vertices.extend(ring2)
    ring3 = [vector + [0] for vector in get_demicube_vertices(6, alt_mode)]
    vertices.extend(ring3)
  return vertices

# bases = get_bases()
bases = get_bases_v2()
# bases = get_bases_rearranged()
# bases = get_e6_bases()
# bases = get_2_21_bases()
# bases = get_1_22_bases()

# v6d = get_cube_vertices(6)
# v6d = get_6_demicube_vertices()
# v6d = get_6_demicube_vertices_alt()
# v6d = get_2_21_vertices()
v6d = get_4_21_vertices_v2()
# v6d = get_4_21_vertices()
# v6d = get_2_31_vertices()
# v6d = get_3_21_vertices()

# v6d = get_1_22_vertices()
v3d = project_to_3d(v6d, bases)

edges = get_edges(v6d)

draw_wireframe(v3d, edges)
# print(len(v6d))
# print(len(edges))

