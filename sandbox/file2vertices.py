"""
--- Currently doesn't work ---
It is supposed to:
- Read in vertices from file
- Find the order of vertices for which the edges don't intersect
- Use this for a BoundaryPolygon object to write DNA inside of

does not work for shapes with concave sections (e.g. a star)
and only for .stl files at the moment
"""

import meshio
import numpy as np
from drawNA.polygons import make_polygon, BoundaryPolygon
import itertools
import math
from shapely.geometry import LineString


def main():
    # read file
    mesh = meshio.read("/home/shanil/drawNA/polygons/STL/hexagon.stl")
    # output vertices
    mesh = mesh.points
    # use only vertices on the z = 0 plane
    mesh = mesh[mesh[:, 2] == 0]

    mesh = np.random.permutation(np.array(mesh))

    permutations = itertools.permutations(range(0, mesh.shape[0]))

    SHAPE = BoundaryPolygon(mesh)
    SHAPE2 = SHAPE
    SHAPE2 = np.vstack((SHAPE.vertices, SHAPE.vertices))

    # make_polygon(mesh)

    print(f"mesh shape is {mesh.shape[0]} \n and SHAPE2 is {SHAPE2}")

    def shapeIsSimple(SHAPE3):
        # make_polygon(SHAPE3)
        return LineString(SHAPE3).is_simple

    perm = next(permutations)

    for j in range(0, math.factorial(mesh.shape[0])):

        SHAPE3 = SHAPE2[perm[0]]

        for i in range(1, mesh.shape[0]):
            SHAPE3 = np.vstack((SHAPE3, SHAPE2[perm[i]]))

        if shapeIsSimple(SHAPE3):  # no intersections
            print(f"Finally \n Updated orientation: {SHAPE3}")
            break

        perm = next(permutations)

    print(f"Permutation: {perm}")

    make_polygon(SHAPE3)


if __name__ == "__main__":
    main()
