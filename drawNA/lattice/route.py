from typing import List

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from ..tools import DNANode, DNAEdge
from .edge import LatticeEdge
from .node import LatticeNode
from ..oxdna import Strand, System
from ..oxdna.strand import generate_helix, POS_BACK, FENE_LENGTH

class LatticeRoute(Strand):
    """
    Strand which follows the path defined by a set of
    LatticeNodes which are joined by LatticeEdges
    """
    def __init__(self, nodes : List[LatticeNode] = []):
        self._nodes = nodes
        self.update_nucleotides()
        super().__init__(self._nucleotides)

    @property
    def nodes(self) -> List[LatticeNode]:
        return self._nodes

    def add_node(self, addition: LatticeNode, index: int = None, update : bool = False):
        """
        Method to add strand(s) to the system

        Parameters:
            addition - accepted as Strand objects or a List of Strands
            index (default = None) - Strand will append to current system,
            otherwise Strand inserted at location given
                
        """

        try:
            assert isinstance(addition, LatticeNode)
        except TypeError:
            raise TypeError(f"addition must be LatticeNode but is {type(addition)}")

        try:
            if index == None:
                self._nodes.append(addition.copy())
            else:
                self._nodes.insert(index, addition.copy())
        except TypeError:
            raise TypeError("Index must an an integer")

        if update:
            self.update_nucleotides()

    def add_nodes(self, node_obj: (list or dict) = None, index: int = None):
        """
        Add multiple strands to the system. Use a list of strands with
        an index indicating where they start, or a dictionary where
        the key is an integer and the value is a strand. 
        
        strand_lists are added in reverse order 
        to preserve original list order

        strand_dicts are added in normal order
        to preserve dictionary keys

        Parameters:
            - node_obj (None) : list or dict of strands
            - index (None) : index to start adding list
        """
        if isinstance(node_obj, list):
            for node in node_obj[::-1]:
                self.add_node(node, index)

        elif isinstance(node_obj, dict):
            for index, node in sorted(node_obj.items()):
                self.add_node(node, index)

        else:
            raise TypeError(
                "add_nodes() requires ONE of a list or dictionary of nodes"
            )

        self.update_nucleotides()
    
    @property
    def edges(self):
        _edges = [LatticeEdge(node, self.nodes[i+1]) for i, node in enumerate(self.nodes[:-1])]
        return _edges

    def update_nucleotides(self):
        self._nucleotides = []
        for i, edge in enumerate(self.edges):
            strand = edge.strand()[0].copy()
            self._nucleotides += strand.nucleotides

    def update_strands(self):
        """
        Generate strands between 2D lattice nodes,
        then generate a single strand representing the route
        """
        strands = []
        edge_0 = self.edges[0]
        strand = generate_helix(n = edge_0.number_of_nt,
                                start_position= edge_0.vertices[0],
                                direction = edge_0.unit_vector,
                                a1 = edge_0.perp_vector)[0]
        strands.append(strand.copy())

        if len(self.edges) < 2:
            return strands            

        for i in range(2,len(self.edges)+1,2):

            last_nuc = strands[-1].nucleotides[-1]
            direction = -last_nuc._a3
            a1 = -last_nuc._a1

            # ensure the backbone position is FENE_LENGTH away from
            # the backbone position of the previous nucleotide
            start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

            # generate strand above that's going in opposite direction
            strand = generate_helix(n=self.edges[i].number_of_nt,
                                    start_position=start,
                                    direction=direction,
                                    a1=a1)[0]
            strands.append(strand)

        return strands

    def get_strand(self):
        nucleotides = []
        _strands = self.update_strands()
        for strand in _strands:
            nucleotides += strand.nucleotides
        _strand = Strand(nucleotides=nucleotides)
        return [_strand, _strands]

    def plot(self, ax : plt.Axes = None, fout : str = None, aspect : int = 5, colour = 'k'):
        nodes = np.array(self.nodes)
        if not ax:
            fig, ax = plt.subplots()
        #plt.grid(True)
        ax.xaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.plot(nodes[:, 0], nodes[:, 1], 'kx', ms = 0.5)
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("No. of strands")
        for edge in self.edges:
            ax.arrow(
                edge.vertices[0][0], 
                edge.vertices[0][1], 
                edge.vector[0], 
                edge.vector[1], 
                width=0.02,
                color = colour,
                length_includes_head=True)
        plt.gca().set_aspect(aspect)
        if fout:
            plt.savefig(fout, dpi=500)
        plt.show()

    def system(self, strands = None, **kwargs) -> System:
        _strand = [self.get_strand()[0]]
        if strands:
            for strand in strands:
                _strand.append(strand)
        _system = System(kwargs.get('box', np.array([50., 50., 50.])))
        _system.add_strands(_strand)

        return _system