from typing import List

import numpy as np

from .oxdna.strand import Strand, generate_helix
from .oxdna.utils import get_rotation_matrix, next_5p, SHIFT_ACROSS, SHIFT_BASE

class DNAObject:
    def __init__(self):
        self.units = 'oxdna'

class DNANode(np.ndarray):
    """
    Abstract class for use with DNAEdge that helps to determine
    angles and vectors needed to generate stable structures.
    """

    def __new__(cls, *args, **kwargs):
        return np.ndarray.__new__(cls, (3)) * 0.0

    def __init__(self, position: np.ndarray, **kwargs):
        self[:] = position[:]
        self._vector_3p = kwargs.get('vector_3p', None)
        self._vector_5p = kwargs.get('vector_5p', None)
        self._pos_3p = kwargs.get('pos_3p', None)
        self._pos_5p = kwargs.get('pos_5p', None)
        self._a1_3p = kwargs.get('a1_3p', None)
        self._a1_5p = kwargs.get('a1_5p', None)
        self._a3_3p = kwargs.get('a3_3p', None)
        self._a3_5p = kwargs.get('a3_5p', None)
        self._angle_3p = 0.0
        self._angle_5p = 0.0

    @property
    def angle(self) -> float:
        """
        Returns the angle between the two vectors
        leaving the Node in radians
        """
        if self._vector_3p is None:
            return None
        if self.vector_5p is None:
            return None

        return np.arccos(np.dot(-self._vector_3p, self._vector_5p))

    @property
    def vector_3p(self) -> np.ndarray:
        """
        Vector entering/leaving the node from the 3' direction

        vector_3p and vector_5p follow on from each other which
        means that they cannot start/finish at the same point
        """
        return self._vector_3p

    @vector_3p.setter
    def vector_3p(self, new_vector: np.ndarray):
        self._vector_3p = new_vector

    @property
    def vector_5p(self) -> np.ndarray:
        """
        Vector leaving/entering the node from the 5' direction

        vector_3p and vector_5p follow on from each other which
        means that they cannot start/finish at the same point
        """
        return self._vector_5p

    @vector_5p.setter
    def vector_5p(self, new_vector: np.ndarray):
        self._vector_5p = new_vector

    @property
    def a1_3p(self):
        return self._a1_3p

    @a1_3p.setter
    def a1_3p(self, new_vector : np.ndarray):
        self._a1_3p = new_vector

    @property
    def a1_5p(self):
        return self._a1_3p

    @a1_5p.setter
    def a1_5p(self, new_vector : np.ndarray):
        self._a1_5p = new_vector

    @property
    def a3_3p(self):
        return self._a1_3p

    @a3_3p.setter
    def a3_3p(self, new_vector : np.ndarray):
        self._a3_5p = new_vector

    @property
    def a3_5p(self):
        return self._a1_3p

    @a3_5p.setter
    def a3_5p(self, new_vector : np.ndarray):
        self._a3_5p = new_vector

    @property
    def pos_3p(self):
        return self._pos_3p

    @pos_3p.setter
    def pos_3p(self, new_vector : np.ndarray):
        self._pos_3p = new_vector

    @property
    def pos_5p(self):
        return self._pos_5p

    @pos_5p.setter
    def pos_5p(self, new_vector : np.ndarray):
        self._pos_5p = new_vector

class DNAEdge:
    """
    Abstract class that allows subclasses to generate oxDNA Strand
    instances along a vector.
    """

    def __init__(
        self, 
        vertex_1: (np.ndarray or DNANode), 
        vertex_2: (np.ndarray or DNANode),
        theta: float = 0.0,
    ):
        self.theta = theta
        if not isinstance(vertex_1, DNANode):
            vertex_1 = DNANode(vertex_1)
        if not isinstance(vertex_2, DNANode):
            vertex_2 = DNANode(vertex_2)
        assert vertex_1 is not vertex_2
        self.vertices = (vertex_1, vertex_2)
        self.vertices[0].vector_5p = self.vector
        self.vertices[0].a3_5p = self.unit_vector
        
        self.vertices[1].vector_3p = -self.vector
        self.vertices[1].a3_3p = self.unit_vector

    def strand(self, sequence: str = None, **kwargs) -> List[Strand]:

        if not sequence:
            # in future version, this will not be so straightforward
            no_of_nucleotides_in_edge = self.nt_length
        else:
            no_of_nucleotides_in_edge = len(sequence)
            if len(sequence) >= self.nt_length:
                print(
                    f"FYI: The Length of `sequence` is longer than the max no. of nucleotides "
                    f"that can be contained within this edge, i.e. {self.nt_length} nucleotides"
                )

        if self.vertices[0].a1_5p:
            a1 = self.vertices[0].a1_5p
        else:
            a1 = self.perp_vector

        start = self.vertices[0] - (SHIFT_ACROSS + SHIFT_BASE) * a1

        strands = generate_helix(
            n=no_of_nucleotides_in_edge,
            sequence=sequence,
            start_position=start,
            a1=a1,
            direction=self.unit_vector,
            **kwargs,
        )
        #self.vertices[0].update_from_nucleotide(strands[0].nucleotides[0], '3p')
        #self.vertices[1].update_from_nucleotide(strands[0].nucleotides[-1], '5p')
        return strands

    def segments(self) -> float:
        """
        There are 2.5 segments in every unit of oxDNA distance.

        Each segment represents a Nucleotide
        """
        return self.length / 2.5

    def node(self, node_3p: DNANode = None, node_5p: DNANode = None) -> DNANode:
        """
        Returns a DNANode for the opposite end of the DNANode provided in 
        parameters
        """
        if not (node_3p or node_5p) or (node_3p and node_5p):
            raise TypeError(
                "Only give one node which is at the" " 3' or 5' end of the Edge"
            )
        if node_3p:
            pass
        elif node_5p:
            pass
        else:
            raise TypeError("Shouldn't get to this point")

    @property
    def summary(self) -> str:
        """Extended equivalent of self.__repr__"""
        string = [""]
        string.append(f"{self.vertices[0]} -> {self.vertices[1]}")
        string.append(f"{self.unit_vector}")
        string.append(f"r={self.length}, nt={self.nt_length}")
        return "\n\t".join(string)

    @property
    def length(self):
        """
        The length of the edge in geometry units
        e.g. for nodes [0,0,0] and [5,0,0] length equals 5
        """
        return np.linalg.norm(self.vector)

    @property
    def nt_length(self):
        """The distance from the first nucleotide to last nucleotide in oxDNA units"""
        return int(self.length * 2.45)

    @property
    def number_of_nt(self):
        """
        The number of nucleotides in an edge object
        e.g. for nodes [0,0,0] and [5,0,0] there would be 6 nucleotides
        """
        return int(np.linalg.norm(self.vector)+1)

    @property
    def vector(self) -> np.ndarray:
        return self.vertices[1] - self.vertices[0]

    @property
    def unit_vector(self) -> np.ndarray:
        return self.vector / (self.vector ** 2).sum() ** 0.5

    @property
    def perp_vector(self) -> np.ndarray:
        """Perpendicular vector which lies in the xy plane"""
        return np.cross(self.unit_vector, np.array([0, 0, 1]))
