from drawNA.lattice import LatticeRoute, LatticeEdge, LatticeNode
from drawNA.oxdna import Strand, System
from drawNA.tools import DNANode, DNAEdge
from drawNA.lattice.utils import find_crossover_locations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
from typing import List

class StapleNode(DNANode):
    def __init__(self, position : np.ndarray):
        super().__init__(position)

class StapleEdge(DNAEdge):
    def __init__(self, vertex_1: StapleNode, vertex_2: StapleNode):
        super().__init__(vertex_1, vertex_2)

class StapleRoute(Strand):
    def __init__(self, scaffold_rows: List[Strand], nodes: List[StapleNode] = []):
        self._nodes = nodes
        self._scaffold_rows = scaffold_rows
        self.update_strand_and_nucleotides()
        super().__init__(self._nucleotides)


    @property
    def nodes(self) -> List[StapleNode]:
        return self._nodes

    @property
    def edges(self):
        _edges = [StapleEdge(node, self.nodes[i+1]) for i, node in enumerate(self.nodes[:-1])]
        return _edges
    
    @property
    def scaffold_rows(self) -> List[Strand]:
        """ Returns scaffold as it's constituent rows"""
        return self._scaffold_rows

    def update_strand_and_nucleotides(self, **kwargs):
        self._nucleotides = []
        for i, edge in enumerate(self.edges):
            if i % 2 == 1:
                continue

            x1 = int(edge.vertices[0][0])
            x2 = int(edge.vertices[1][0])
            row = int(edge.vertices[0][1])
            if len(self.scaffold_rows) - 1 == row:
                break
            nucleotides = []
            scaffold_row = self.scaffold_rows[row]
            
            if x1 > x2:
                x1, x2 = x2, x1
                for nucleotide in scaffold_row.nucleotides[x1:x2+1]:
                    nucleotides.append(nucleotide.make_across())
                # switch the order of the nucleotides back again
                nucleotides = nucleotides[::-1]
            else:
                for nucleotide in scaffold_row.nucleotides[::-1][x1:x2+1]:
                    nucleotides.append(nucleotide.make_across())
 
            for nucleotide in nucleotides:
                self._nucleotides.append(nucleotide)



class StapleCollection:
    """ Probably need to edit this such that
    it becomes easy to join staples on the same row together
    not sure, how to implement that yet """
    def __init__(self, strands: List[Strand] = []):
        self._strands = strands
    
    @property
    def staples(self) -> List[Strand]:
        return self._strands

    @property
    def n_staples(self) -> int:
        return len(self._strands)

    def add_staples(self, staple_strand: Strand):
        self._strands.append(staple_strand)

    def plot_nodes(self, strand: Strand, ax, colour = 'r', width = 0.02, **kwargs):
        nodes = np.array(strand.nodes)
        #plt.grid(True)
        ax.plot(nodes[:, 0], nodes[:, 1], 'bx', ms = 0.5)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("No. of strands")
        for edge in strand.edges:
            ax.arrow(
                edge.vertices[0][0], 
                edge.vertices[0][1], 
                edge.vector[0], 
                edge.vector[1], 
                width=width,
                color = colour,
                length_includes_head=True, **kwargs)
    

    def plot(self, route: LatticeRoute = None, fout: str = None):
        fig, ax = plt.subplots()
        if route:
            self.plot_nodes(strand = route, ax = ax, colour = 'k', width = 0.1, alpha = 0.)
        for staple in self.staples:
            colour = np.random.rand(3,)
            self.plot_nodes(strand = staple, ax = ax, colour = colour)
        
        plt.gca().set_aspect(5)
        if fout:
            plt.savefig(fout, dpi=500)

        plt.show()



    

class ScaffoldRows:
    def __init__(self,route: LatticeRoute, staple_widths = [5, 16, 26]):
        self.route = route
        self._staple_widths = staple_widths

    @property
    def sizes(self) -> list:
        """ Returns length of all the rows"""
        _row_sizes = []
        for i in range(0,len(self.route.edges),2):
            _row_sizes.append(self.route.edges[i].nt_length + 1)
        return _row_sizes

    @property
    def bounds(self) -> list:
        """ Returns x coord of start and end point """
        _bounds = []
        for i in range(0,len(self.route.edges),2):
            start = self.route.edges[i].vertices[0][0]
            end = self.route.edges[i].vertices[1][0]
            _bounds.append([start,end])
        return _bounds

    @property
    def start_side(self) -> list:
        """ Returns left or right """
        _start_side = []
        for i in range(0,len(route.nodes), 2):
            side = "left" if route.nodes[i][0]-route.nodes[i+1][0] < 0 else "right"
            _start_side.append(side)
        return _start_side

    def n_staples(self, staple_width = None):
        """ Returns no. of staples per row """
        if not staple_width: # default to 16
            staple_width = self._staple_widths[1]
        return [int(i/staple_width) for i in self.sizes]
    
    def unpaired_bases(self, staple_width = None):
        """ Returns no. of unpaired bases per row """
        if not staple_width: # default to 16
            staple_width = self._staple_widths[1]
        return [int(i % staple_width) for i in self.sizes]

    @property
    def info(self, staple_width = None):
        """
        Returns information for each row in the Lattice Route
        Specifically its size and start/end point
        And the number of staples 
        """
        return pd.DataFrame({
                "size": self.sizes,
                "bounds": self.bounds,
                "start side": self.start_side,
                "staples (5)": self.n_staples(self._staple_widths[0]),
                "staples (16)": self.n_staples(self._staple_widths[1]),
                "unpaired bases (5)": self.unpaired_bases(self._staple_widths[0]),
                "unpaired bases (16)": self.unpaired_bases(self._staple_widths[1]),
            })
        
def side_staples(route : LatticeRoute, staple_width = 16):
    scaffold = ScaffoldRows(route)
    staples = []
    scaffold_rows = route.get_strand()[1]
    for row, info in scaffold.info.iterrows():
        x2 = info["bounds"][0]
        if info["start side"] == "left":
            x1 = x2 + staple_width
        elif info["start side"] == "right":
            x1 = x2 - staple_width
        
        staple_nodes = [
            StapleNode([x1, row, 0.0]),
            StapleNode([x2, row, 0.0]),
            StapleNode([x2, row+1, 0.0]),
            StapleNode([x1, row+1, 0.0])]

        staples.append(StapleRoute(scaffold_rows, staple_nodes))

    return StapleCollection(staples)
    

def test_StapleRoute(route: LatticeRoute):
    x1 = 0
    x2 = 15
    row = 0
    scaffold_rows = route.get_strand()[1]
    staple_nodes = [
        StapleNode([x1, row, 0.0]),
        StapleNode([x2, row, 0.0]),
        StapleNode([x2, row+1, 0.0]),
        StapleNode([x1, row+1, 0.0])]
    
    return StapleRoute(scaffold_rows, staple_nodes)


if __name__ == "__main__":
    half_turn_indices   = [4, 15, 25, 36, 46, 56, 67, 77, 88, 98, 109]
    staple_lengths      = [9, 31, 51, 73]

    route_vertices = [
        [0., 0.,0.],
        [56.,0.,0.],
        [56.,1.,0.],
        [0., 1.,0.],
        [0., 2.,0.],
        [56.,2.,0.],
        [56., 3.,0.],
        [0., 3.,0.],
        [0., 4.,0.],
        [56., 4.,0.],
        [56.,5.,0.],
        [0., 5.,0.],
        [0., 6.,0.],
        [56.,6.,0.],
        [56., 7.,0.],
        [0., 7.,0.],
        [0., 8.,0.],
        [88., 8.,0.],
        [88., 9.,0.],
        [0., 9.,0.],
        [0., 10.,0.],
        [88., 10.,0.],
        [88., 11.,0.],
        [0., 11.,0.],
        [0., 12.,0.],
        [56., 12.,0.],
        [56., 13.,0.],
        [0., 13.,0.],
        [0., 14.,0.],
        [56., 14.,0.],
        [56, 15, 0],
        [0,15,0]
    ]

    nodes = [LatticeNode(np.array(i)) for i in route_vertices]
    route = LatticeRoute(nodes)
    fig, ax = plt.subplots()
    route.plot(ax = ax)


    # test, scaf = test_StapleRoute(route)
    # system = System(np.array([50,50,50]))
    # system.add_strands(scaf)
    # system.write_oxDNA("scaffold")
    collection = side_staples(route)
    system = route.system(collection.staples)
    system.write_oxDNA("lol")
    