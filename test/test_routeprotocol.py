from drawNA.oxdna import System
from drawNA.lattice import LatticeRoute
from drawNA.polygons import BoundaryPolygon

import numpy as np
import matplotlib.pyplot as plt
from os import path

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])

DIST_SQUARE = 2.60

def square():
    square = np.array([[0, 0, 0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    polygon = BoundaryPolygon(square * 20)
    return polygon

def square_route():
    polygon = square()
    lattice = polygon.dna_snake(
        grid_size=(DIST_SQUARE, 5.44)
    )
    route = lattice.route()
    return route

def test_square():
    route = square_route()
    system = System(np.array([50.0, 50.0, 50.0]))
    system.add_strand(route)
    system.write_oxDNA(root=ROOT)

if __name__ == "__main__":
    polygon = square()
    print(polygon.vertices)
    route = square_route()
    fig, ax = plt.subplots()
    polygon.plot2D(ax=ax)
    route.plot(ax=ax)
    plt.show()
    test_square()
