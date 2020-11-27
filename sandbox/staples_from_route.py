from drawNA.lattice import LatticeRoute, LatticeEdge, LatticeNode
from drawNA.oxdna import Strand, System
from drawNA.tools import DNANode, DNAEdge
from drawNA.lattice.utils import find_crossover_locations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
from typing import List
import plotly.express as px

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

class StapleBaseClass:
    """
    Before the staples are generated, we need to investigate the scaffold path (lattice route)
    and process this into information which we can use to design staples across the structure.
    """
    def __init__(self, scaffold: LatticeRoute):
        self.route = scaffold
        self.scaffold_rows = self.route.get_strand()[1]
        self.n_rows = len(self.scaffold_rows)

        # So that the values do not need to be re-computed everytime their are called
        self._row_seq = self.get_row_sequence()
        self._start_side = self.get_start_side()
        self._row_size = self.get_sizes()
        self._bounds = self.construct_bounds()

        self.lattice = self.make_lattice()
        print(f"StapleBaseClass generated.")

    @property
    def row_size(self) -> List[int]:
        """ Returns length of all the rows"""
        return self._row_size
    
    @property
    def row_seq(self) -> List[str]:
        "Returns sequence of each row, given in the 3' >> 5' direction"
        return self._row_seq

    @property
    def start_side(self) -> List[str]:
        """ 
        Returns which side is the 3' for each row, i.e.:
        3' >> 5' = left   or  5' << 3' = right
        """
        return np.array(self._start_side)

    @property
    def bounds(self) -> list:
        """ Returns x coord of start (3') and end (5') of each row """
        return np.array(self._bounds)

    def start(self) -> list:
        """ Returns x coord of start (3' side) """
        return self._bounds[:,0]

    def end(self) -> list:
        """ Returns x coord of end (5' side) """
        return self._bounds[:,1]

    @property
    def unpaired_nt(self) -> List[int]:
        unpaired_nt = []
        for row in range(self.n_rows):
            unpaired_nt.append(row, np.count_nonzero(hello.lattice[row,:,0]==0))   
        return np.array(unpaired_nt)
        
    @property
    def info(self):
        """
        Returns information for each row in the Lattice Route
        Specifically its size and start/end point
        And the number of staples 
        """
        return pd.DataFrame({
                "size": self.row_size,
                "bounds": list(self.bounds),
                "start side": self.start_side,
                "no. of unpaired nt": self.unpaired_nt
            })
    
    

    def make_lattice(self) -> np.ndarray:
        """
        This generates a 3D numpy array (n x rows x 2)
        star
        """
        # Create a correctly sized lattice filled with 0's
        width = int(self.bounds.max())
        lattice = np.zeros((self.n_rows,width,2))

        # sort so index 0 is smallest (left bound) and 1 is therefore the right bound
        bounds = np.sort(self.bounds)

        # For every row, find left and right bounds
        # if left != 0, make all points before left bound = None
        # if right != width of array, make all points after left bound = None
        for row in range(self.n_rows):
            if bounds[row,0] != 0: # self.bounds.min()
                lattice[row, 0:bounds[row,0]] = [None, None]

            if bounds[row,1] != width:
                lattice[row, bounds[row,1]::] = [None, None]
        print("Initial Staple Lattice Formed")
        return lattice

    def plot_lattice(self):
        fig = px.imshow(self.lattice[:,:,0], aspect=4)
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='White')
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='White')
        fig.show()
        return None

    def get_sizes(self) -> list:
        print("Retrieving Row Sizes")
        """ Returns length of all the rows"""
        _row_sizes = []
        for i in range(0,len(self.route.edges),2):
            _row_sizes.append(self.route.edges[i].nt_length)

            # Assertion each edge should equal length of sequence in scaffold rows
            assert _row_sizes[int(i/2)] == len(self.row_seq[int(i/2)])
        return _row_sizes

    def get_row_sequence(self):
        print("Retrieving Row Sequences")
        "Returns sequence (for the 3' >> 5' direction) for each row"
        row_seq = []
        for i in range(0,self.n_rows):
            row_seq.append(self.scaffold_rows[i].sequence)
        return row_seq

    def get_start_side(self) -> list:
        print("Retrieving Start Side")
        """ 
        Returns which side is the 3' for each row, i.e.:
        3' >> 5' = left   or  5' << 3' = right
        """
        _start_side = []
        route = self.route
        for i in range(0,len(route.nodes), 2):
            side = "left" if route.nodes[i][0]-route.nodes[i+1][0] < 0 else "right"
            _start_side.append(side)
        return _start_side

    def construct_bounds(self) -> List[List[int]]: #as np.ndarray
        print("Retrieving Bounds")
        """ Returns x-coord of start (3' side) and end (5' side) nodes """
        _bounds = []

        ## Find x-coords
        x_pos = 0.0 # Initialise to 0x
        for row in range(self.n_rows):
            if self.start_side[row] == "left":
                x_3p = x_pos
                x_5p = x_pos + self.row_size[row]

            elif self.start_side[row] == "right":
                x_3p = x_pos 
                x_5p = x_pos - self.row_size[row]

            else:
                print("Something has gone very wrong, this should not happen")
            
            _bounds.append([x_3p, x_5p])

            x_pos = x_5p

        ## Normalise bounds to start at 0
        _bounds = np.array(_bounds, dtype=int)
        _bounds -= _bounds.min()
        assert _bounds.min() == 0
        return _bounds

        
        # return np.array(_bounds)
        # for i in range(0,len(self.route.edges),2):
        #     start = self.route.edges[i].vertices[0][0]
        #     end = self.route.edges[i].vertices[1][0]
        #     _bounds.append([start,end])
        
        # # Normalise Bounds to start at 0
        # _bounds = np.array(_bounds, dtype=int)
        # _bounds -= _bounds.min()
        # assert _bounds.min() == 0
        # return _bounds



class StaplingAlgorithm1(StapleBaseClass):
    def __init__(self, scaffold, option2=None, option3=None):
        super.__init__(self, scaffold)

    

        
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
        [-40., 1.,0.],
        [-40., 2.,0.],
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
    nodes = [LatticeNode(np.array(i)) for i in np.array(route_vertices)/[5,1,1]]
    route = LatticeRoute(nodes)
    # rows = ScaffoldRows(route)
    
#     fig, ax = plt.subplots()
#     route.plot(ax = ax)

    hello = StapleBaseClass(route)
    display(hello.info)
    hello.plot_lattice()


    
    



# Migrate To New File "test_staplebaseclass.py"

# def test_makelattice():
#     nodes = [
#         [0.,0.,0.],
#         [60.,0.,0.],
#         [60.,1.,0.],
#         [40.,1.,0.],
#         [40.,2.,0.],
#         [80.,2.,0.],
#         [80.,3.,0.],
#         [0.,3.,0.],
#         [0.,4.,0.],
#         [60.,4.,0.]
#     ]
#     route_nodes = [LatticeNode(np.array(i)) for i in nodes]
#     route = LatticeRoute(route_nodes)
#     hello = [print(i) for i in route.edges[i].summary]
#     staple_base = StapleBaseClass(route)
#     base_lattice = staple_base.make_lattice()
#     return route, base_lattice, hello, staple_base

# if __name__ == "__main__":
#     route, lattice, summary, staple_base = test_makelattice()

    # test, scaf = test_StapleRoute(route)
    # system = System(np.array([50,50,50]))
    # system.add_strands(scaf)
    # system.write_oxDNA("scaffold")
    # collection = side_staples(route)
    # system = route.system(collection.staples)
    # system.write_oxDNA("lol")



    ### Final Workflow should be something like
    # collections = system.generate_staples(options) # here system should only contain a single scaffold strand
        # - collections consists of 


