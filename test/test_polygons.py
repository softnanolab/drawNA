#!/usr/bin/env python
"""
pytest module for drawNA.polygons
"""
from os import path

import numpy as np

from drawNA.polygons import Edge, BoundaryPolygon
from drawNA.oxdna import Strand

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])


def test_Edge():
    edge = Edge(np.array([0.0, 1.0, 0.0]), np.array([0.0, -5.5, 0.0]))
    assert edge.length == 6.5
    assert edge.kind == "boundary"
    assert np.linalg.norm(edge.vector - np.array([0, -6.5, 0.0])) == 0.0
    assert np.linalg.norm(edge.unit_vector - np.array([0, -1, 0])) == 0.0

    helix = edge.strand(sequence="AAAGGG")
    print(helix)
    assert helix[0].sequence == "AAAGGG"


def test_BoundaryPolygon():
    square = BoundaryPolygon(np.array([[1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0],]))

    assert square.n_edges == len(square.edges)
    assert square.x.sum() == 2
    assert square.y.sum() == 2
    assert square.z.sum() == 0
    assert isinstance(square.edges[0], Edge)

    # square.plot2D(show=False)
    square.write_STL(f"{ROOT}/square.stl")
    square.write_PLY(f"{ROOT}/square.ply")


if __name__ == "__main__":
    test_Edge()
    test_BoundaryPolygon()
