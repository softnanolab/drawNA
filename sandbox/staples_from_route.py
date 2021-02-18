from drawNA.lattice import LatticeRoute, LatticeEdge, LatticeNode
from drawNA.oxdna import Strand, System
from drawNA.tools import DNANode, DNAEdge
from drawNA.lattice.utils import find_crossover_locations
from drawNA.polygons import BoundaryPolygon

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd

from copy import deepcopy
from typing import List

import itertools as it

from softnanotools.logger import Logger

logger = Logger("Staples_from_route")
logger.level = 10

class StapleBaseClass:
    """
    Before the staples are generated, we need to investigate the scaffold path (lattice route)
    and process this into information which we can use to design staples across the structure.
    """
    def __init__(self, scaffold: LatticeRoute, crossover_threshold: float):
        print(f"StapleBaseClass: Initialising StapleBaseClass object {self}")
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
        self.staple_ID = 1
        print(f"StapleBaseClass generated.\n-----")
    
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
        3' to 5' = left   or  5' from 3' = right
        """
        return np.array(self._start_side)

    @property
    def bounds(self) -> list:
        """ Returns x coord of start (3') and end (5') of each row """
        return np.array(self._bounds)

    @property
    def unpaired_nt(self) -> List[int]:
        """ Returns the number of unpaired nucleotides on each row of the scaffold"""
        unpaired_nt = []
        for row in range(self.n_rows):
            # count number of sites with a value of zero (as opposed to the name of the function)
            unpaired_nt.append(np.count_nonzero(self.lattice[row,:,0]==0))

        return np.array(unpaired_nt)
    
    @property
    def total_unpaired_nt(self) -> List[int]:
        """ Returns the total number of unpaired nucleotides in the configuration"""
        return np.sum(self.unpaired_nt)

    def start(self) -> list:
        """ Returns x coord of start (3' side) """
        return self._bounds[:,0]

    def end(self) -> list:
        """ Returns x coord of end (5' side) """
        return self._bounds[:,1]


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
    
    def make_lattice(self, threshold) -> np.ndarray:
        """
        This generates a 3D numpy array (max no. of nucleotides x no. of scaffold rows x 4 layers)

        Firstly, there are 4 layers (it is a 3D array [width, height, depth=4])
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
            
            This is important for converting the staple array to a container of staples
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
        Layer B stores all crossover positions along each scaffold helix
        i.e. 1 (crossover up) -1 (crossover down)
        
        By filtering through the a1 vector of all the nucleotides, those furthest from the
        central helical axes in the XY plane (i think) are deemed the "crossover" positions.

        The a1 vector seems to lie in the range of -1.09 < a1 > 1.09 [roughly]

        Parameters:
            lattice (np.array): contains 0 (empty lattice sites) and None (invalid lattice sites)
            threshold (int): y distance from central helical axis 

        Yields a lattice where:
            1   (crossover up)
            -1  (crossover down)
            0   (non crossover positions) [unchanged]
            None (non scaffold positions) [unchanged]
        
        """

        ## PSEUDO CODE
        # repeat for all the strands in the scaffold
        # for a strand, find the a1 vector of each nucleotide
            # if this vector is greater than the +ve threshold value then return the x position and "U"
            # if this vector is smaller than the -ve threshold value then return the x position and "D"
            # update the row of the lattice representing that strand with U and D at specified x positions from the starting point of that strand in the lattice
        layerB = lattice
        
        for row_N, strand in enumerate(self.scaffold_rows):
            # Find nucleotide distance from central axis of the helix
            distance_from_central_axis = []
            distance_from_central_axis = [nuc._a1[1] for nuc in strand.nucleotides]
            # print(f"Row: {row_N}, Distance from Central: {len(distance_from_central_axis)}, Number of zeros: {np.count_nonzero(layerB[row_N,:]== 0)}")
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

    def plot_lattice(self,layer = 0, aspect = 6, show_staple_nodes = False):
        lattice_2D = self.lattice[:,:, layer]
        fig, ax = plt.subplots(figsize=(15,15))
        # fig.set_facecolor("grey")
        vals = np.linspace(0,1,256)
        np.random.shuffle(vals)
        cmap = plt.cm.colors.ListedColormap(plt.cm.viridis(vals))
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

        plt.gca().set_aspect(aspect)
        plt.gca().invert_yaxis()
        ax.set_title("Crossover positions")
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("Scaffold row strands")
        plt.show()
        return

    def fill_lattice_with_single_domains(self):
        """Adds single domained staples to all un-paired nucleotide sections in the lattice 
        
        T0D0 Future: 
            - add the ability to staple sections of each row, i.e. if the staples aren't at the edges, there would be 2-3 staples on one row
            - only do this where there is a sufficient number of nucleotides
        """

        lattice = self.lattice[:,:,0]
        info = self.info_dataframe

        for row, lattice_row in enumerate(lattice):
            unstapled_nt_indices = np.where(lattice_row == 0)[0]

            if info['no. of unpaired nt'][row] <= 2:
                continue
            elif info['start side'][row] == 'right':
                staple_node_S = [unstapled_nt_indices[0], row, 'S']
                staple_node_T = [unstapled_nt_indices[-1], row, 'T']
            else:
                staple_node_S = [unstapled_nt_indices[-1], row, 'S']
                staple_node_T = [unstapled_nt_indices[0], row, 'T']
            
            staple = [staple_node_S, staple_node_T]
            self.write_staples_to_array([staple])

    def write_staples_to_array(self, staples: List[List]):
        """ Records the position of staple(s) nt and their start/stop points in the lattice

        Parameters:
            staples (list) - a list of staples, where each `staple` 
                has an even number of `nodes` (list items)
                
                Each `node` is a 3 item list, in the form [<x index>,<row index>,<staple type>]
                These specify a particular nucleotide on a particular row in the scaffold. The 
                staple type is either S, C or T; standing for start, crossover or terminus, 
                respectively. The start will be the 3p side and 5p side will be the terminus.
        
        Example input:
            staple_1 = [[0,0,'S'],[10,0,'T']] # Single domain staple
            staple_2 = [[10,1,'S'],[0,1,'C'],[0,2,'C'],[10,2,'T']] # Double domain staple
            write_staples_to_array([staple_1, staple_2])

        """
        assert type(staples) == list
        assert len(staples) <= 2
        staple_type_to_int = {'S': 10, 'C': 20, 'T':30}
        for staple in staples:
            # Ensure staple is valid
            if not any(x[0] is None for x in staple):

                # assign type of staple ('S' 'T' or 'C') to layer C
                for node in staple:
                    [x, row, staple_type] = node 
                    self.lattice[row, x, 2] = staple_type_to_int[staple_type]
                
                # assign staple_ID to lattice A(0) for all staple nt's
                for i in range(0, len(staple),2):
                    node_1 = staple[i]
                    node_2 = staple[i+1]
                    domain = np.sort([ node_1[0], node_2[0] ])
                    print(f"Writing staple #{self.staple_ID} with domain: {domain} on row {row}")
                    assert node_1[1] == node_2[1]
                    row = node_1[1]
                    self.lattice[row, domain[0]:domain[1]+1, 0] = self.staple_ID
                
                self.staple_ID += 1
            # print(f"Written staples into array. ID: {self.staple_ID}")
            else:
                pass

    def generate_container(self):
        """ Converts `self.lattice` to staples using `self.scaffold_obj` returning `StapleContainer` object

        Generates staples in the following steps:
        - Cycle through staples by staple_ID.
        - First, find the indices of all the points where staple_ID exists
        - Filter indices through layer 2, where values equal 10, 20 or 30 (start, crossover, terminus)
        - Append the scaffold nucleotide index to each index (i.e. row, x-index)
        - Order indices based on row
        - Order indices such that S->C and C->T
        - Order pairs of indices such that S->C...C->T
        - Retrieves corresponding scaffold indices between the staple nodes (i.e. S, C, T)
        - Generates staple and adds to list
        - Converts list of staples to a container
        """
        layer0 = self.lattice[:,:,0] # Staple_IDs
        layer2 = self.lattice[:,:,2] # Special nucleotide: Start (10), Crossovers (20) or End (30). Otherwise (0).
        layer3 = self.lattice[:,:,3] # Scaffold Indices

        # Dictionary converter, float to staple type: Start, Crossover, End
        flt_2_st = {0.0: 0.0, 10.0: 'S', 20.0: 'C', 30.0: 'T'}
        staples = []

        for staple_ID in np.arange(1, self.staple_ID):

            # nested list of indices
            # e.g. [ [0,0], [0,1], [0,2]...[0,21], [1,0], [1,1], [1,2]...[1,21] ]
            staple_indices = np.argwhere(layer0 == staple_ID) 

            # finds staple_type from layer2 and filters list to nt's of special type (S, C or T)
            # e.g. [[ 0,0,'T' ], [ 0,21,'C' ], [ 1,21,'C' ], [ 1,0,'S' ], [1, 21, 'C']]
            staple_indices_and_type = [ [int(row), int(x_index), flt_2_st[layer2[row, x_index]] ] 
                                        for [row, x_index] in staple_indices if layer2[row, x_index] != 0]

            # finds scaffold indices from layer3 and appends to each item in the list
            # e.g. [[0, 0, 'T', 0.0], [0, 21, 'C', 21.0], [1, 0, 'S', 365.0], [1, 21, 'C', 344.0]]
            staple_nodes_info = [ [ row, x_index, staple_type, int(layer3[row, x_index]) ]
                                    for [row, x_index, staple_type] in staple_indices_and_type]

            # e.g. [['1', '0', 'S', '365'], ['1', '21', 'C', '344'], ['0', '21', 'C', '21'], ['0', '0', 'T', '0']]
            ordered_staple_nodes = self._order_staple_indices(staple_nodes_info)
            
            # e.g. [365, 344, 21, 0]
            scaffold_node_indices = [int(info[3]) for info in ordered_staple_nodes]

            # e.g. [[365, 344], [21, 0]]
            no_of_domains = int(len(scaffold_node_indices)/2)
            scaffold_domains = np.array(scaffold_node_indices).reshape(no_of_domains,2)

            # [365, 363, ...345, 344, 21, 19, ... 1, 0]
            all_scaffold_indices = self._make_indices(scaffold_domains)

            staple_nuc = []
            for index in all_scaffold_indices:
                staple_nuc.append(self.scaffold_obj.nucleotides[index].make_across())
            
            staple = Strand(nucleotides=staple_nuc)
            staples.append(staple)
        
        return StapleContainer(staple_strands = staples, staple_base_class = self)
    
    @staticmethod
    def _order_staple_indices(arr: List[list]) -> list:
        """ Orders the staple indices for generating staples from the scaffold

        The order will be organised such that the starting nucleotide (S) will be first between
        each pair of nucleotides. And conversely, the terminating nucleotide (T) will be last
        between each pair of nucleotides. 
        
        The nucleotides which are adjacent to a nucleotide on a different scaffold row 
        (i.e. the corners of a staple) are the crossovers (C).

        Arg:
            arr 

        Returns (list):
            Such that S->C and C->T (or S->T for single domained staples)
            AND SC is followed by CT (if it is not a single domained staple)

        """
        # Ensure S is followed by C and C is followed by T
        for i in np.arange(0,len(arr),2):
            if arr[i+1][2] == 'S' or arr[i][2] == 'T':
                arr[i], arr[i+1] = arr[i+1], arr[i]

        # Ensure S is the first point of the first domain
        if arr[0][2] != 'S':
            arr = np.roll(arr,2,axis=0).tolist()

        # Ensure T is the last point
        assert arr[-1][2] == 'T', f"Likely due to error in staple generation where two staples intersect: {arr}"

        return arr 

    @staticmethod
    def _make_indices(domain_pairs: List[list]) -> list:
        """Makes a single list of integers depicting scaffold indices, in the correct order.

        Arg:
            domain_pairs (list(list)) - e.g. [ [0,10], [40,30] ]

        Returns:
            Using the example above it would be [0,1...10,40,39......31,30]

            More importantly a range of {10 to 0} would be returned as 10,9,8...1,0
        """
        scaffold_indices = []
        for [num1, num2] in domain_pairs:            
            scaffold_indices += [i for i in np.linspace(num1, num2, num=abs(num1-num2)+1, dtype = int)]
        return scaffold_indices

class StaplingAlgorithm1(StapleBaseClass):
    """
    Max number of domains: 2
    """
    def __init__(self, scaffold, domain_size=14, crossover_threshold=0.956):
        super().__init__(scaffold, crossover_threshold)
        self.domain_size = domain_size
        print(f"StapleAlgorithm1: Generating side staples {self}")
        self.generate_side_staples()
        print(f"StapleAlgorithm1: Generating inside staples of single domain only {self}")
        self.generate_inside_staples()
    
    def generate_side_staples(self):
        """ Generates staple nodes for staples at the left and right edges of a scaffold

        T0D0: Currently the workflow doesn't allow for the find_crossovers function
        to check whether a single domained staple is being placed on a site where a
        pre-existing staple sits...
        
        Yields:
            Up to two sets of staple nodes for every two rows
        """
        df = self.info_dataframe

        for row1_idx in range(0, self.n_rows, 2):
            
            # The number of rows is odd, the last row will be single domained
            if row1_idx == (self.n_rows - 1):
                single_domain = True
                # print(f"single domain: {row1_idx, self.n_rows}")
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
                        staple_direction = 'up')

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

                # register staples in the array
                self.write_staples_to_array([staple_nodes_left, staple_nodes_right])
                
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
                # register staples in the array
                self.write_staples_to_array([staple_nodes_left])
                
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

                # register staples in the array
                self.write_staples_to_array([staple_nodes_right])
            

    def generate_inside_staples(self):
        """ Generates single domain staples for all remaining unpaired nt's """
        self.fill_lattice_with_single_domains()

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

            The site types in layer B of the lattice are:
                0: no crossover
                1: upwards pointing scaffold crossover
                -1: downwards pointing scaffold crossover

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
                # THen we find the next suitable crossover point (moving right in our search)
                crossover_indices = np.where(self.lattice[row_index, :, 1] == target_scaffold_direction)[0]
                x_index  = crossover_indices[crossover_indices > x_index][0]
                
                scaf_site_direction = self.lattice[row_index, x_index, 1]
                assert scaf_site_direction == -staple_direction, "Searching for staple didn't work"

                # Domain size is set to be one nt smaller than the location where there is a crossover
                if self.lattice[row_index, x_index, 0] != 0:

                    # Find location of pre-existing crossover(s) [20: 'C']
                    # Find the one which is at the index < x_index
                    existing_crossover_indices = np.where(self.lattice[row_index, :, 2] == 20)[0]
                    x_index  = existing_crossover_indices[existing_crossover_indices <= x_index][-1] - 1

            elif self.lattice[row_index, x_index, 2] == 20: # checking for a crossover already existing at this spot
                x_index -= 1
            
            # Ensure distance between bound and x_index is positive 
            if x_index - bound > 0:
                return int(x_index)
            # Staples with None will not be written in the `write_staples_to_array` function
            else:
                return None
        
        else: #direction is right to left
        
            x_index = bound - self.domain_size

            scaf_site_direction = self.lattice[row_index, x_index, 1]

            if scaf_site_direction != target_scaffold_direction:
                
                # Finding index of all the sites on the row where the staple
                # will have the correct direction (opposite to scaffold direction)
                # THen we find the next suitable crossover point (moving left in our search)
                crossover_indices = np.where(self.lattice[row_index, :, 1] == target_scaffold_direction)[0]
                x_index  = crossover_indices[crossover_indices < x_index][-1]

                scaf_site_direction = self.lattice[row_index, x_index, 1]
                assert scaf_site_direction == -staple_direction, "Searching for staple didn't work"

                # Domain size is set to be one nt smaller than the location where there is a crossover
                if self.lattice[row_index, x_index, 0] != 0:
                    
                    # Find location of pre-existing crossover(s) [20: 'C']
                    # Find the one which is at the index > x_index
                    existing_crossover_indices = np.where(self.lattice[row_index, :, 2] == 20)[0]
                    x_index  = existing_crossover_indices[existing_crossover_indices >= x_index][0] + 1

            elif self.lattice[row_index, x_index, 2] == 20: # checking for a crossover already existing at this spot
                x_index += 1

            # Ensure distance between bound and x_index is negative 
            if x_index - bound < 0:
                return int(x_index)
            # Staples containing None will not be written in the `write_staples_to_array` function
            else:
                return None

class StaplingAlgorithm2(StapleBaseClass):
    """
    This algorithm is specifically for rectangular shapes. i.e. the edges must be straight

    It goes about putting 4 staples on each row. 2 from the left, 2 from the right (shifted up one row).

    Steps:
    1. Calculate approximate index of positions of 5 potential start/crossover/end points
        for staples added to the scaffold
        - L  (First nt) 
        - Q1 (nt at approx the 1/4 mark)
        - Q2 (nt at approx the 2/4 mark)
        - Q3 (nt at approx the 3/4 mark)
        - R  (Last nucleotide)
        
    2. Do an oscillatory search around those points for the coordinates 
        corresponding to a crossover in the scaffold
    3. 
    """
    def __init__(self, scaffold, crossover_threshold=0.956):
        super().__init__(scaffold, crossover_threshold)
        self.verify_route_is_applicable()

        self.L = self.Q1 = self.Q2 = self.Q3 = self.R = 0
        self.find_approx_crossover_indices()
        self.find_closest_crossover_indices()
        
        self.generate_staples()
        print(self.info_dataframe)

    def verify_route_is_applicable(self):
        """Checks to see left and right edges of scaffold are aligned"""
        _3p_idx = self.info_dataframe['3p']
        _5p_idx = self.info_dataframe['5p']
        min, max = _3p_idx[0], _5p_idx[0]
        all_idx = list(_3p_idx) + list(_5p_idx)
        print("StapleAlgorithm2: Verifying inputted route is supported")
        if not all(idx == min or idx == max for idx in all_idx):
            logger.kill(f'Scaffold route is unsupported with this algorithm, edges must be aligned')

    def find_approx_crossover_indices(self):
        """ Calculates approximate indices of 5 potential crossover points """
        pass

    def find_closest_crossover_indices(self):
        """ Ensures 5 potential crossover points all occur where  """
        pass
    
    def generate_staples(self):
        """ Generates staples in two cycles: 1) left half, 2) right half """
        pass

class StapleContainer:
    """ Stores staple and scaffold strand objects with the ability to modify them.

    T0D0: add methods to modify staples -> could be added through subclassing
    """
    def __init__(self, 
                staple_strands: List[Strand] = [], 
                staple_base_class: StapleBaseClass = None):

        self._staples = staple_strands
        self._baseclass = staple_base_class
        self._scaffold = self._baseclass.route
    
    @property
    def staples(self) -> List[Strand]:
        """List of staples as `Strand` objects"""
        return self._staples

    @property
    def scaffold(self) -> LatticeRoute:
        """The Scaffold `Strand` object as a `LatticeRoute` object
        
            Note: `LatticeRoute` is a subclass of `Strand`
        """
        return self._scaffold
    
    @property
    def base_class(self) -> StapleBaseClass:
        return self._baseclass

    def system(self, **kwargs) -> System:
        _system = System(kwargs.get('box', np.array([50., 50., 50.])))
        _system.add_strands(self.staples)
        _system.add_strand(self.scaffold)
        return _system
    
    @property
    def n_staples(self) -> int:
        return len(self._staples)

    def add_staples(self, staple_strand: Strand):
        self._staples.append(staple_strand)

    def plot(self):
        self._baseclass.plot_lattice()


## Copied from protocols/lattice-route/DNA_snake.py
def generate(polygon_vertices: np.ndarray, title: str = "Generate()") -> LatticeRoute:
    print('Running Generate Function')
    print(f'Creating polygon from {len(polygon_vertices)} vertices')
    polygon = BoundaryPolygon(polygon_vertices)
    print("Constructing scaffolding lattice")
    # print(f"{title}: ...constructing scaffolding lattice...")
    lattice = polygon.dna_snake(start_side="left", grid_size = [0.34, 2])
    # print(f"{title}: ...calculating route.")
    print("Calculating Route")
    route = lattice.route()
    print("Plotting Route")
    route.plot()
    return route

## Protocol functions
def staple_1_and_write_to_file(route: LatticeRoute, name_of_file: str, domain_size = 15):
    print("Class StaplingAlgorithm1: Generating side staples")
    staple_1 = StaplingAlgorithm1(route, domain_size = domain_size)
    # staple_1.fill_lattice_with_single_domains()

    print("staple_1_and_write_to_file(): Adding staples to a container...")
    container_1 = staple_1.generate_container()
    
    print("staple_1_and_write_to_file(): Adding staples to an oxDNA system...")
    system_1 = container_1.system()
    
    # print("staple_1_and_write_to_file(): Writing `.top` and `.conf` files")
    # system_1.write_oxDNA(prefix = name_of_file)
    return system_1, container_1

def plot_staples(staples: StapleContainer):
    if type(staples) == StapleContainer:
        staples = staples.base_class

    print("Plotting Staples with start, crossover and terminal points labelled")
    staples.plot_lattice(layer=0, show_staple_nodes=True)
    print("Plotting Staples with scaffold crossover points labelled")
    staples.plot_lattice(layer=0, show_staple_nodes=False)

def main():
    hourglass = np.array([[0.,0.,0.],[4,6.,0.],[0,12,0],[12,12,0],[8,6.,0.],[12,0.,0.]])
    stacked_I = np.array([
    [0.,0.,0.],[3.,0.,0.],[3.,1.,0.],[2.,1.,0.], [2.,2.,0.],[3.,2.,0.],
    [3.,3.,0.],[2.,3.,0.],[2.,4.,0.],[3.,4.,0.],[3.,5.,0.],[0.,5.,0.],[0.,4.,0.],[1.,4.,0.],
    [1.,3.,0.],[0.,3.,0.], [0.,2.,0.],[1.,2.,0.],[1.,1.,0.],[0.,1.,0.]
    ])
    square = np.array([[0,0,0],[10,0,0],[10,9,0],[0,9,0]])
    triangle = np.array([[0,0,0],[5,10,0],[10,0,0]])
    trapREV = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]])

    route = generate(square*2)
    # staple, container = staple_1_and_write_to_file(route, "square25", domain_size=25)
    # plot_staples(container)
    staple_2 = StaplingAlgorithm2(route)
    return staple_2

if __name__ == "__main__":
    staple_2 = main()
    
    # route = generate(stacked_I*17)
    # system, container = staple_and_write_to_file(route, "stacked_I")
    # plot_staples(container)

    # route = generate(hourglass*10)
    # system, container = staple_and_write_to_file(route, "hourglass")
    # plot_staples(container)

    # route = generate(square*10)
    # system, container = staple_and_write_to_file(route, "square")
    # plot_staples(container)

    # route = generate(trapREV*17)
    # system, container = staple_and_write_to_file(route, "trapezium")
    # plot_staples(container)

    # route = generate(square*2)

    # container = staple_and_write_to_file(route, "square25", domain_size=25)
    # plot_staples(container)
    # container = staple_and_write_to_file(route, "square15", domain_size=15)
    # plot_staples(container)
    # container = staple_and_write_to_file(route, "square35", domain_size=35)
    # plot_staples(container)
    # container = staple_and_write_to_file(route, "square45", domain_size=45)
    
    # plot_staples(container)
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
    # container = side_staples(route)
    # system = route.system(container.staples)
    # system.write_oxDNA("lol") 



    ### Final Workflow should be something like
    # containers = system.generate_staples(options) # here system should only contain a single scaffold strand
        # - containers consists of 


# import numpy as np
# import matplotlib.pyplot as plt

# staple_left = [[0,0,'S'],[8,0,'C'],[8,1,'C'],[0,1,'T']]
# staple_right = [[19,1,'S'],[12,1,'C'],[12,0,'C'],[19,0,'T']]
# staple_type_to_int = {'S': 10, 'C': 20, 'T':30}
# staple_ID = 2

# def write_staples_to_array(staples: list, lattice, staple_ID = 1):
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
# grid = write_staples_to_array([staple_left, staple_right], grid)
# plt.imshow(grid[:,:,0])
# plt.gca().set_title("Purple=unpaired, otherwise, stapled")
# plt.show()
# plt.imshow(grid[:,:,2])
# plt.show()