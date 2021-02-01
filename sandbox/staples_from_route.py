from drawNA.lattice import LatticeRoute, LatticeEdge, LatticeNode
from drawNA.oxdna import Strand, System
from drawNA.tools import DNANode, DNAEdge
from drawNA.lattice.utils import find_crossover_locations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
import pandas as pd
from typing import List
import plotly.express as px
from drawNA.polygons import BoundaryPolygon
import math
from copy import deepcopy

import itertools as it

class StapleNode(DNANode):
    def __init__(self, position : np.ndarray):
        super().__init__(position)

class StapleEdge(DNAEdge):
    def __init__(self, vertex_1: StapleNode, vertex_2: StapleNode):
        super().__init__(vertex_1, vertex_2)

class StapleRoute(Strand):
    """T0D0 - Re-write this function"""
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
        """List of staples as `Strand` objects"""
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
    def __init__(self, scaffold: LatticeRoute, crossover_threshold: float):
        self.route = scaffold
        self.scaffold_obj, self.scaffold_rows = self.route.get_strand()
        self.n_rows = len(self.scaffold_rows)

        # Define here so not to re-compute values everytime they are accessed
        self._row_seq = self.get_row_sequence()
        self._start_side = self.get_start_side()
        self._row_size = self.get_number_of_nt()
        self._bounds = self.construct_bounds()

        self.lattice_width = int(self.bounds.max()+1) # +1 to account for counting starting at 0
        self.lattice = self.make_lattice(threshold = crossover_threshold)
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
        """ Returns the number of unpaired nucleotides on each row of the scaffold"""
        unpaired_nt = []
        for row in range(self.n_rows):
            # count number of sites with a value of zero (as opposed to the name of the function lol)
            unpaired_nt.append(np.count_nonzero(self.lattice[row,:,0]==0))
        return np.array(unpaired_nt)
        
    @property
    def info_dataframe(self):
        """Returns scaffold geometry information.

        Returns information for each row in the Lattice Route.
        Specifically its size, start/end point, which side the 3' begins
            and the number of unpaired nt        

        T0D0: The number of staples 
        """
        return pd.DataFrame({
                "size": self.row_size,
                "3p": self.bounds[:,0],
                "5p": self.bounds[:,1],
                "start side": self.start_side,
                "no. of unpaired nt": self.unpaired_nt
            })
    
    def help_lattice(self):
        """How is this constructed?
        
        Firstly, there are 3 layers (it is a 3D array [width, height, depth=2])
            Where cells have the value 'None', means that no scaffold lives here
        
        Layer A: Contains the Staple ID (if exists) otherwise 0 for empty lattice site
            
            StapleID: i.e. 1,2,3,4,5,6,... 
                a way to separate different staples
            
            Unhybridised lattice site: 0
                a way to determine whether a staple has been generated for 
                the particular site on the scaffold

            Invalid lattice site: NoneType
        
        Layer B: Details all scaffold crossover locations
            i.e. 1 (crossover up) -1 (crossover down)

        Layer C: Contains info about start, end and turns of the staple

            Empty space: 0 
            --- the same as layer 1
            
            Start of a staple: 10
            --- 3' side of staple
            
            Crossovers: 20
            
            End of a staple: 30
            --- 5' side of the staple


        Layer D: Contains index information of the scaffold.
            
            This is important for converting the staple array to a collection of staples
        """
    
    def make_lattice(self, threshold) -> np.ndarray:
        """
        This generates a 3D numpy array (n x rows x 2)
        star
        """
        # Create a correctly sized lattice filled with 0's
        lattice = np.zeros((self.n_rows,self.lattice_width))

        # Sort bounds
        # left bound (index 0) is smallest
        # therefore, right bound (index 1) is biggest
        bounds = np.sort(self.bounds)

        # For every row, find left and right bounds
        # all points before left bound = None
        # all points after right bound = None
        for row in range(self.n_rows):
            
            if bounds[row,0] != 0: # self.bounds.min()
                lattice[row, 0:bounds[row,0]] = None

            if bounds[row,1] != self.lattice_width - 1: # to account for bounds counting from 0
                lattice[row, bounds[row,1]+1::] = None

        layerA = deepcopy(lattice) # 0 (empty space) and 1,2,3... (staple ID) or None (invalid lattice site)
        layerB = self.make_layerB(deepcopy(lattice), threshold) # crossover loc (1: 'up', -1: 'down')
        layerC = deepcopy(lattice) # Start, end and turns of the staple
        layerD = self.make_layerD(deepcopy(lattice))

        lattice_3D = np.stack([layerA,layerB,layerC, layerD],axis=2) 
        return lattice_3D

    def make_layerB(self, lattice, threshold):
        """
        Layer B contains all crossover positions along each scaffold helix
        i.e. 1 (crossover up) -1 (crossover down)
        
        Parameters:
            lattice (np.array): contains 0 (empty lattice sites) and None (invalid lattice sites)
            threshold (int): y distance from central helical axis 
        
        """

        ## PSEUDO CODE
        # for a strand, find the a1 vector of each nucleotide
        # if this vector is greater than the +ve threshold value then return the x position and "U"
        # if this vector is smaller than the -ve threshold value then return the x position and "D"
        # update the row of the lattice representing that strand with U and D at specified x positions from the starting point of that strand in the lattice
        # repeat for all the strands in the scaffold
        layerB = lattice
        
        for row_N, strand in enumerate(self.scaffold_rows):
            # Find nucleotide distance from central axis of the helix
            distance_from_central_axis = []
            distance_from_central_axis = [nuc._a1[1] for nuc in strand.nucleotides]
            print(f"Row: {row_N}, Distance from Central: {len(distance_from_central_axis)}, Number of zeros: {np.count_nonzero(layerB[row_N,:]== 0)}")
            assert len(distance_from_central_axis) == np.count_nonzero(layerB[row_N,:]== 0)

            # Find indices where crossover is possible
            arr = np.array(distance_from_central_axis)
            cross_up_idx = np.where(arr > threshold)[0]
            cross_down_idx = np.where(arr < -threshold)[0]


            if self.start_side[row_N] == "left":
                # Start side is found as the first 0 from the left side
                start_idx = np.where(lattice[row_N]==0)[0][0]
                # Shift crossover indices relative to start index
                cross_up_idx += start_idx
                cross_down_idx += start_idx

            else: # start_side = right
                # Start side is found as the first 0 from the right side
                start_idx = np.where(lattice[row_N]==0)[0][-1]
                # Shift crossover indices relative to start index
                # i.e. we must count these backwards, i.e. [0, 10] -> [56, 46]
                cross_up_idx = start_idx - np.array(cross_up_idx)
                cross_down_idx = start_idx - np.array(cross_down_idx)

            # Modify layerB to represent crossover locations
            # cross_up_idx and cross_down_idx are arrays representing the crossover indices
            # 1  -> crossover up
            # -1 -> crossover down
            np.put(layerB[row_N], cross_up_idx,    1)
            np.put(layerB[row_N], cross_down_idx, -1)
        
        return layerB

    def make_layerD(self, lattice) -> np.ndarray:
        """
        Layer D contains information of what the nucleotide
        indices are for a certain scaffold route.

        This will be very important when we want to convert
        the staple array into staple objects.
        """

        start = 0
    
        for row_number, size in enumerate(self.row_size):
            assert(np.count_nonzero(lattice[row_number]==0)) == size 
            
            # Index of the final nucleotide on `row_number`
            end = start + size
            array_of_indices = np.arange(start,end)
            
            # Flip indices for rows with 3' on the right
            if self.start_side[row_number] == 'right':
                array_of_indices = array_of_indices[::-1]

            # Replace lattice sites (0) with the index array
            np.place(
                arr = lattice[row_number],
                mask = lattice[row_number] == 0,
                vals = array_of_indices
            )
            
            start = end

        return lattice

    def plot_lattice(self,layer = 1, aspect = 6, show_staple_nodes = False):
        lattice_2D = self.lattice[:,:, layer]
        fig, ax = plt.subplots(figsize=(15,15))
        # fig.set_facecolor("grey")
        cmap = cm.get_cmap("coolwarm")
        cmap.set_bad(color='black')
        plt.imshow(lattice_2D, cmap=cmap,aspect = 2)

        if show_staple_nodes == False:
            lattice_crossovers = self.lattice[:,:, 1]
            for i in range(lattice_crossovers.shape[0]):
                for j in range(lattice_crossovers.shape[1]):
                    if lattice_crossovers[i, j] == -1:
                        text = "↓"#"▼""D"
                    elif lattice_crossovers[i, j] == 1:
                        text = "↑"#"▲""U"
                    else:
                        text = ""
                    ax.text(j, i, text,ha="center", va="center", color="w")
        else:
            lattice_crossovers = self.lattice[:,:, 2]
            for i in range(lattice_crossovers.shape[0]):
                for j in range(lattice_crossovers.shape[1]):
                    if lattice_crossovers[i, j] == 10:
                        text = "S"#"▼""D"
                    elif lattice_crossovers[i, j] == 20:
                        text = "C"#"▲""U"
                    elif lattice_crossovers[i, j] == 30:
                        text = "T"#"▲""U"
                    else:
                        text = ""
                    ax.text(j, i, text,ha="center", va="center", color="w")
        ## For layerD
        # lattice_crossovers = self.lattice[:,:, 1]
        # for i in range(lattice_crossovers.shape[0]):
        #     for j in range(lattice_crossovers.shape[1]):
        #         text = lattice_crossovers[i, j]
        #         ax.text(j, i, text,ha="center", va="center", color="w", fontsize=5, rotation = 90)

        plt.gca().set_aspect(aspect)
        plt.gca().invert_yaxis()
        ax.set_title("Crossover positions")
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("Scaffold row strands")
        plt.show()
        return None

    def get_number_of_nt(self) -> list:
        print("Retrieving Row Sizes")
        """ Returns length of all the rows"""
        _row_sizes = []
        for i in range(0,len(self.route.edges),2):
            _row_sizes.append(self.route.edges[i].number_of_nt)

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
                x_5p = x_pos + self.row_size[row] - 1 

            elif self.start_side[row] == "right":
                x_3p = x_pos
                x_5p = x_pos - self.row_size[row] + 1

            else:
                print("Something has gone very wrong, this should not happen")
            
            _bounds.append([x_3p, x_5p])

            x_pos = x_5p

        ## Normalise bounds to start at 0
        _bounds = np.array(_bounds, dtype=int)
        _bounds -= _bounds.min()
        assert _bounds.min() == 0
        return _bounds

    def lattice_to_staples(self):
        """ Converts `self.lattice` to staples using `self.scaffold_obj`

        Generates staple 
        
        """
        pass

class StaplingAlgorithm1(StapleBaseClass):
    """
    Max number of domains: 2
    """
    def __init__(self, scaffold, domain_size=14, crossover_threshold=0.956):
        super().__init__(scaffold, crossover_threshold)
        self.domain_size = domain_size
        self.staple_ID = 1
        print(self.info_dataframe)
        self.generate_side_staples()
    
    def generate_side_staples(self):

        df = self.info_dataframe

        for row1_idx in range(0, self.n_rows, 2):
            
            # The number of rows is odd, the last row will be single domained
            if row1_idx == (self.n_rows - 1):
                single_domain = True
                print(f"single domain: {row1_idx, self.n_rows}")
            else:
                single_domain = False

            """Make a single domained staple"""
            if single_domain:
                row1 = df.loc[row1_idx]

                if row1['start side']=='left': # 3p -------> 5p
                    boundary_node_left = row1['3p']
                    crossover_node_left = self.find_crossover(
                        bound = boundary_node_left,
                        row_index = row1_idx,
                        direction = 'left_to_right',
                        staple_direction = 'down')

                    boundary_node_right = row1['5p']
                    crossover_node_right = self.find_crossover(
                        bound = boundary_node_right,
                        row_index = row1_idx,
                        direction = 'right_to_left',
                        staple_direction = 'up')

                    ## Staple directions
                    # Scaffold 3p --------------------> 5p
                    # Staples  SN2<--SN1       SN2<--SN1
                    staple_node_1L = [crossover_node_left, row1_idx, 'S']
                    staple_node_2L = [boundary_node_left, row1_idx, 'T']

                    staple_node_1R = [boundary_node_right, row1_idx, 'S']
                    staple_node_2R = [crossover_node_right, row1_idx, 'T']

                else: # 5p <------- 3ps
                    boundary_node_left = row1['5p']
                    crossover_node_left = self.find_crossover(
                        bound = boundary_node_left,
                        row_index = row1_idx,
                        direction = 'left_to_right',
                        staple_direction = 'down')

                    boundary_node_right = row1['3p']
                    crossover_node_right = self.find_crossover(
                        bound = boundary_node_right,
                        row_index = row1_idx,
                        direction = 'right_to_left',
                        staple_direction = 'down')

                    ## Staple directions
                    # Scaffold 5p <-------------------- 3p
                    # Staples  SN1-->SN2       SN1-->SN2
                    staple_node_1L = [boundary_node_left, row1_idx, 'S']
                    staple_node_2L = [crossover_node_left, row1_idx, 'T']

                    staple_node_1R = [crossover_node_right, row1_idx, 'S']
                    staple_node_2R = [boundary_node_right, row1_idx, 'T']
                
                # Assign staples
                staple_nodes_left = [staple_node_1L, staple_node_2L]
                staple_nodes_right = [staple_node_1R, staple_node_2R]
                
            else:
                ### ASSIGN ROW 1 AND ROW 2 DATAFRAMES
                row2_idx = row1_idx + 1
                row1 = df.loc[row1_idx]
                row2 = df.loc[row2_idx]

                ### FINDING CROSSOVER POSITIONS OF THE STAPLE
                """ ADD STAPLES TO LEFT SIDE"""

                # Row 1: 3p ---------------------> 5p
                if row1['start side'] == 'left':
                    boundary_node_row_1_L = row1['3p']
                    boundary_node_row_2_L = row2['5p']

                    #  Left edge is aligned both row1 and row2 have the same left bound
                    #  Or Left edge is MISaligned: row2 bound is further left than row1
                    if (row1['3p'] == row2['5p']) or (row1['3p'] > row2['5p']):
                        crossover_node_left = self.find_crossover(
                            bound = boundary_node_row_1_L,
                            row_index = row1_idx, 
                            direction = 'left_to_right',
                            staple_direction = 'up')
                        
                    # Left edge is MISaligned: row1 is further left than row2
                    elif row2['5p'] > row1['3p']:
                        crossover_node_left = self.find_crossover(
                            bound = boundary_node_row_2_L,
                            row_index = row2_idx, 
                            direction = 'left_to_right',
                            staple_direction = 'down')
                    
                    ## Schematic of staple 
                    # Row 2   (3p)  SN1----->SN2 
                    # Row 1   (5p)  SN4<-----SN3
                    # Note: crossover between SN2 and SN3
                    staple_node_1L = [boundary_node_row_2_L, row2_idx, 'S']
                    staple_node_2L = [crossover_node_left, row2_idx, 'C']
                    staple_node_3L = [crossover_node_left, row1_idx, 'C']
                    staple_node_4L = [boundary_node_row_1_L, row1_idx, 'T'] 

                # Row 1: 5p <--------------------- 3p
                else:
                    boundary_node_row_1_L = row1['5p']
                    boundary_node_row_2_L = row2['3p']

                    # left edge, by definition, must be aligned
                    if row1['5p'] == row2['3p']:
                        crossover_node_left = self.find_crossover(
                            bound = boundary_node_row_1_L,
                            row_index = row1_idx, 
                            direction = 'left_to_right',
                            staple_direction = 'up')
                    else:
                        print("This should not be possible")

                    ## Schematic of staple 
                    # Row 2   (5p)  SN4<-----SN3 
                    # Row 1   (3p)  SN1----->SN2
                    # Note: crossover between SN2 and SN3
                    staple_node_1L = [boundary_node_row_1_L, row1_idx, 'S'] 
                    staple_node_2L = [crossover_node_left, row1_idx, 'C']
                    staple_node_3L = [crossover_node_left, row2_idx, 'C']
                    staple_node_4L = [boundary_node_row_2_L, row2_idx, 'T'] 

                staple_nodes_left = [staple_node_1L, staple_node_2L, staple_node_3L, staple_node_4L]

                """ADD STAPLES TO THE RIGHT SIDE"""

                # Row 1: 5p <---------------------- 3p
                if row1['start side'] == 'right':
                    boundary_node_row_1_R = row1['5p']
                    boundary_node_row_2_R = row2['3p']

                    # Row 2: 3p ------------> 5p             3p ------------> 5p
                    # Row 1: 5p <------------ 3p     OR      5p <------ 3p
                    if (row1['3p'] == row2['5p']) or (row1['3p'] < row2['5p']):
                        crossover_node_right = self.find_crossover(
                            bound = boundary_node_row_1_R,
                            row_index = row1_idx, 
                            direction = 'right_to_left',
                            staple_direction = 'up')
                        
                    # Row 2: 3p ------> 5p
                    # Row 1: 5p <----------- 3p
                    elif row2['5p'] < row1['3p']:
                        crossover_node_right = self.find_crossover(
                            bound = boundary_node_row_2_R,
                            row_index = row2_idx, 
                            direction = 'right_to_left',
                            staple_direction = 'down')
                    
                    ## Schematic of staple 
                    # Row 2   (3p)  SN2<-----SN1 
                    # Row 1   (5p)  SN3----->SN4
                    # Note: crossover between SN2 and SN3
                    staple_node_1R = [boundary_node_row_2_R, row2_idx, 'S']
                    staple_node_2R = [crossover_node_right, row2_idx, 'C']
                    staple_node_3R = [crossover_node_right, row1_idx, 'C']
                    staple_node_4R = [boundary_node_row_1_R, row1_idx, 'T'] 

                # Row 1: 3p ---------------------> 5p
                else:
                    boundary_node_row_1_R = row1['5p']
                    boundary_node_row_2_R = row2['3p']

                    # Row 2: 3p <------------ 3p  
                    # Row 1: 5p ------------> 5p   i.e. right side is a crossover
                    if row1['5p'] == row2['3p']:
                        crossover_node_right = self.find_crossover(
                            bound = boundary_node_row_1_R,
                            row_index = row1_idx, 
                            direction = 'right_to_left',
                            staple_direction = 'up')
                    else:
                        print("This should not be possible")

                    ## Schematic of staple 
                    # Row 2   (5p)  SN3----->SN4 
                    # Row 1   (3p)  SN2------SN1
                    # Note: crossover between SN2 and SN3
                    staple_node_1R = [boundary_node_row_1_R, row1_idx, 'S'] 
                    staple_node_2R = [crossover_node_right, row1_idx, 'C']
                    staple_node_3R = [crossover_node_right, row2_idx, 'C']
                    staple_node_4R = [boundary_node_row_2_R, row2_idx, 'T'] 

                staple_nodes_right = [staple_node_1R, staple_node_2R, staple_node_3R, staple_node_4R]

            """REGISTER NODES"""
            self.write_nodes_to_array([staple_nodes_left, staple_nodes_right])
            

    def generate_inside_staples(self, fixed_domain_length : int = 15):
        # cycle through the rows in chunks (number of domains)
        pass



    def find_crossover(self,
            bound: int,
            row_index: int,
            direction: str,
            staple_direction: str
            ) -> int:
        """Returns x index for the crossover point of a staple.

        T0D0: add a self.staple_exists() function and call it before returning value

        Args:
            bound (int) - the start or end x-index of a staple

            row_index (int) - which row of the scaffold we are using
                to find a crossover point

            direction (str) - from which side, relative to the given `bound`,
                a crossover should be found
                i.e. "left_to_right" or "right_to_left"
            
            staple_direction (str) - which direction the staple will
                crossover to the next helix
                i.e. 'up' or 'down'

        Returns:
            The x-index of a crossover position for a staple.

            This is specifically chosen to the locations where
            the nucleotide is at a peak position from the central
            helical axes.

            Where the crossover required is up, the scaffold must
            have a peak going down. This corresponds to a value of -1
            in `layerB` of the `lattice`.

        """

        assert direction in ['left_to_right', 'right_to_left']
        assert staple_direction in ['up', 'down']
        
        # If staple crossover is intended to go up, scaffold crossover must be going down (-1)
        crossover_type = {'up' : 1, 'down' : -1}

        staple_direction = crossover_type[staple_direction]
        target_scaffold_direction = -staple_direction

        # Recall: layerB of the lattice is the one denoting positions of crossovers

        ## FIND INDEX
        if direction == 'left_to_right':
            x_index = bound + self.domain_size

            # Site type: 0 - no crossover, -1/1 - crossover position
            scaf_site_direction = self.lattice[row_index, x_index, 1]

            if scaf_site_direction != target_scaffold_direction:
                
                # Finding index of all the sites on the row, where the staple
                # will have the correct direction (opposite to scaffold direction)
                crossover_indices = np.where(self.lattice[row_index, :, 1] == target_scaffold_direction)[0]

                # Domain size is one nucleotide smaller
                if np.any(crossover_indices == x_index - 1):
                    x_index -= 1
                # Domain size is one or more nucleotides greater
                else:
                    x_index  = crossover_indices[crossover_indices > x_index][0]

                scaf_site_direction = self.lattice[row_index, x_index, 1]
            
            # Site type should be -1 or 1. opposite of staple direction
            assert scaf_site_direction == -staple_direction, "Finding Crossovers Failed"
        
        else: #direction is right to left
        
            x_index = bound - self.domain_size

            # Site type: 0 - no crossover, -1/1 - crossover position
            scaf_site_direction = self.lattice[row_index, x_index, 1]

            if scaf_site_direction != target_scaffold_direction:
                
                # Finding index of all the sites on the row where the staple
                # will have the correct direction (opposite to scaffold direction)
                crossover_indices = np.where(self.lattice[row_index, :, 1] == target_scaffold_direction)[0]

               
                if np.any(crossover_indices == x_index + 1):  # Domain size is one nucleotide smaller
                    x_index += 1
                else: # Domain size is one or more nucleotides greater
                    # Finds next suitable crossover point
                    x_index  = crossover_indices[crossover_indices < x_index][-1] # Most right index of crossover points left of the current x_index

                scaf_site_direction = self.lattice[row_index,x_index, 1]

            # Site type should be -1 or 1. opposite of staple direction
            assert scaf_site_direction == -staple_direction, "Finding Crossovers Failed"

        return int(x_index)

    def write_nodes_to_array(self, staples: list):
        assert type(staples) == list
        assert len(staples) <= 2
        staple_type_to_int = {'S': 10, 'C': 20, 'T':30}
        for staple in staples:
            for node in staple:
                # assign type of staple ('S' 'T' or 'C') to lattice C(2)
                [x, row, staple_type] = node 
                # print(f"Node: {node}")
                self.lattice[row, x, 2] = staple_type_to_int[staple_type]
            
            for i in range(0, len(staple),2):
                node_1 = staple[i]
                node_2 = staple[i+1]
                # assign staple_ID to lattice A(0) for all staple nt's
                domain = np.sort([ node_1[0], node_2[0] ])
                print(f"Domain: {domain}")
                assert node_1[1] == node_2[1]
                row = node_1[1]
                self.lattice[row, domain[0]:domain[1]+1, 0] = self.staple_ID
            
            self.staple_ID += 1
            print(f"Written staples into array. ID: {self.staple_ID}")


## Copied from protocols/lattice-route/DNA_snake.py
def generate(polygon_vertices: np.ndarray, title: str = "Generate()"):
    print(f"{title}: Making polygon...")
    polygon = BoundaryPolygon(polygon_vertices)
    print(f"{title}: ...constructing scaffolding lattice...")
    lattice = polygon.dna_snake(straightening_factor=5, start_side="left", grid_size = [0.34, 2])
    print(f"{title}: ...calculating route.")
    route = lattice.route()
    return route

if __name__ == "__main__":
    hourglass = np.array([[0.,0.,0.],[4,6.,0.],[0,12,0],[12,12,0],[8,6.,0.],[12,0.,0.]])
    stacked_I = np.array([
    [0.,0.,0.],[3.,0.,0.],[3.,1.,0.],[2.,1.,0.], [2.,2.,0.],[3.,2.,0.],
    [3.,3.,0.],[2.,3.,0.],[2.,4.,0.],[3.,4.,0.],[3.,5.,0.],[0.,5.,0.],[0.,4.,0.],[1.,4.,0.],
    [1.,3.,0.],[0.,3.,0.], [0.,2.,0.],[1.,2.,0.],[1.,1.,0.],[0.,1.,0.]
    ])
    square = np.array([[0,0,0],[10,0,0],[10,10,0],[0,10,0]])
    route = generate(hourglass*5)
    route.plot()
    staple_1 = StaplingAlgorithm1(route, crossover_threshold=0.956, domain_size=17)
    staple_1.plot_lattice(layer=3)
    display(staple_1.info_dataframe)
    staple_1.plot_lattice(layer=0, show_staple_nodes=True)


    # half_turn_indices   = [4, 15, 25, 36, 46, 56, 67, 77, 88, 98, 109]
    # staple_lengths      = [9, 31, 51, 73]

    # route_vertices = [
    #     [0., 0.,0.],
    #     [56.,0.,0.],
    #     [56.,1.,0.],
    #     [-40., 1.,0.],
    #     [-40., 2.,0.],
    #     [56.,2.,0.],
    #     [56., 3.,0.],
    #     [0., 3.,0.],
    #     [0., 4.,0.],
    #     [56., 4.,0.],
    #     [56.,5.,0.],
    #     [0., 5.,0.],
    #     [0., 6.,0.],
    #     [56.,6.,0.],
    #     [56., 7.,0.],
    #     [0., 7.,0.],
    #     [0., 8.,0.],
    #     [88., 8.,0.],
    #     [88., 9.,0.],
    #     [0., 9.,0.],
    #     [0., 10.,0.],
    #     [88., 10.,0.],
    #     [88., 11.,0.],
    #     [0., 11.,0.],
    #     [0., 12.,0.],
    #     [56., 12.,0.],
    #     [56., 13.,0.],
    #     [0., 13.,0.],
    #     [0., 14.,0.],
    #     [56., 14.,0.],
    #     [56, 15, 0],
    #     [0,15,0]
    # ]
    # nodes = [LatticeNode(np.array(i)) for i in np.array(route_vertices)/[5,1,1]]
    # route = LatticeRoute(nodes)
    # rows = ScaffoldRows(route)

    
#     fig, ax = plt.subplots()
#     route.plot(ax = ax)

    # hello = StapleBaseClass(route)
    # display(hello.info)
    # hello.plot_lattice()
    
    
    



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


# import numpy as np
# import matplotlib.pyplot as plt

# staple_left = [[0,0,'S'],[8,0,'C'],[8,1,'C'],[0,1,'T']]
# staple_right = [[19,1,'S'],[12,1,'C'],[12,0,'C'],[19,0,'T']]
# staple_type_to_int = {'S': 10, 'C': 20, 'T':30}
# staple_ID = 2

# def write_nodes_to_array(staples: list, lattice, staple_ID = 1):
#     assert type(staples) == list
#     assert len(staples) <= 2
#     for staple in staples:
#         print(staple)
#         for node in staple:
#             # assign type of staple ('S' 'T' or 'C') to lattice C
#             [x, row, staple_type] = node 
#             lattice[row, x, 2] = staple_type_to_int[staple_type]
        
#         for i in range(0, len(staple),2):
#             node_1 = staple[i]
#             node_2 = staple[i+1]

#             # assign staple_ID to layer A for all staple nt's
#             domain = np.sort([ node_1[0], node_2[0] ])
#             print(f"Domain: {domain}")
#             row = node_1[1]
#             print(domain[0], domain[1], row)
#             lattice[row, domain[0]:domain[1]+1, 0] = staple_ID
        
#         staple_ID += 1
        
#     return lattice

# grid = np.zeros((5,20,3))
# plt.imshow(grid[:,:,0])
# plt.show()
# grid = write_nodes_to_array([staple_left, staple_right], grid)
# plt.imshow(grid[:,:,0])
# plt.gca().set_title("Purple=unpaired, otherwise, stapled")
# plt.show()
# plt.imshow(grid[:,:,2])
# plt.show()