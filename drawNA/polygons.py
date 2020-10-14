import re
from typing import List

import numpy as np
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from .oxdna import Nucleotide, System, Strand, generate_helix
from .tools import DNAEdge, DNANode
from .lattice import Lattice, DNASnake

EDGE_TYPE = {0: "boundary", 1: "scaffold", 2: "staple"}

class Vertex(DNANode):
    """
    Geometric point in space.
    """
    def __init__(self, position : np.ndarray):
        super().__init__(position)

class Edge(DNAEdge):
    """ 
    An Edge is a combination of two vertices and stores geometric information about this connection 
    """

    def __init__(self, vertex_1: list, vertex_2: list, edge_kind: int = 0):
        # don't use type because that's a protected function in Python
        self.kind = EDGE_TYPE[edge_kind]

        if vertex_1[2] != vertex_2[2]:
            raise ValueError("Edge must lie in the xy plane (z values must be equal)")
        elif (vertex_1 == vertex_2).all():
            raise ValueError("Given verticies must be different")

        # Initialise rest of self from parent class DNAEdge
        super().__init__(Vertex(vertex_1), Vertex(vertex_2))

class BoundaryPolygon:
    """
    A BoundaryPolygon is a cyclic combination of 3 or more edges, which
    result in a closed Polygon

    The first vertex on the starting edge = the last vertex on the final edge
    """

    def __init__(self, vertices: np.ndarray):
        """ 4 or more vertex objects given, each value on a given row """

        assert vertices.shape[1] == 3 or vertices.shape[1] == 2
        self.vertices = vertices

    @property
    def n_edges(self) -> np.array:
        return self.vertices.shape[0]

    @property
    def edges(self) -> List[Edge]:
        """
        Iterate over all vertices and return a list of Edge instances
        """
        edges = []
        for i in range(self.n_edges):
            edges.append(Edge(self.vertices[i - 1], self.vertices[i], 0))
        return edges

    @property
    def x(self) -> np.ndarray:
        return self.vertices[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.vertices[:, 1]

    @property
    def z(self) -> np.ndarray:
        try:
            return self.vertices[:, 2]
        except IndexError:
            raise IndexError(f"Trying to access Polygon.z but {self} is only 2D!")

    def __repr__(self) -> str:
        return f"<Polygon{self.vertices.shape[1]}D Vertices[{self.vertices.shape[0]}]>"

    def __str__(self) -> str:
        return "This polygon has {} edges".format(self.n_edges)

    def write_STL(self, fout: str):
        triangulation = tri.Triangulation(self.vertices[:, 0], self.vertices[:, 1])
        triangles = triangulation.get_masked_triangles()
        cells = [("triangle", triangles)]
        meshio.write_points_cells(fout, self.vertices, cells)

    def write_PLY(self, fout: str, comments: List[str] = []):
        """
        Writes the BoundaryPolygon to a ASCII PLY File.
        """

        # begin header
        output_string = "PLY\nformat ascii 1.0\n"
        for comment in comments:
            output_string += f"comment {comment}\n"
        output_string += f"element vertex {self.vertices.shape[0]}\n"
        output_string += "".join([f"property float {i}\n" for i in ["x", "y", "z"]])
        output_string += (
            "element face 1\nproperty list uchar int vertex_index\nend_header\n"
        )
        # end of header

        # vertex position list
        # regex substitutes out '[' and ']' characters and the replace
        # removes the space at the beginning of each new line
        output_string += re.sub(r"\[|\]", "", self.vertices.__str__()).replace(
            "\n ", "\n"
        )

        # only one face in face list e.g. for a square 4 0 1 2 3
        output_string += f"\n{self.vertices.shape[1]} "
        output_string += " ".join([str(i) for i in range(self.vertices.shape[1])])
        output_string += "\n"

        # finished generating file so now writing
        with open(fout, "w") as f:
            f.write(output_string)

    def plot2D(
        self, ax: plt.Axes = None, fout: str = None, **kwargs
    ):
        """
        Assumes that the shape is 2D and lies on the z=0 plane.

        """
        if not ax:
            fig, ax = plt.subplots()
        ax.plot(list(self.x) + [self.x[0]], list(self.y) + [self.y[0]], "k-", **kwargs)
        ax.set_aspect("equal", "datalim")
        if fout:
            fig.savefig(fout)
        if not ax:
            fig.show()

    def lattice(self, **lattice_kwargs) -> Lattice:
        return Lattice(self.vertices, **lattice_kwargs)

    def dna_snake(self, **lattice_kwargs) -> Lattice:
        return DNASnake(polygon_vertices = self.vertices, **lattice_kwargs)
