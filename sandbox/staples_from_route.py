from drawNA.lattice import LatticeRoute
from drawNA.oxdna import Nucleotide, Strand, System
from drawNA.oxdna import nucleotide
from drawNA.polygons import BoundaryPolygon

from drawNA.oxdna.nucleotide import POS_BASE, SHIFT_BASE 

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib.ticker import MultipleLocator
import pandas as pd

from copy import deepcopy
from typing import List, Tuple

import random

from softnanotools.logger import Logger
from os import path

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

        # Troubleshooting when trialing different scaffold twists...
        # for row_N in range(len(self.scaffold_rows)):
        #     if row_N%2 == 0: #even
                
        #         if layerB[row_N][-1] == 0.0:

        #             layerB[row_N][-2] = 0.0
        #             layerB[row_N][-1] = 1.0
        #     else:
                
        #         if layerB[row_N][-1] == 0.0:
        #             layerB[row_N][-2] = 0.0
        #             layerB[row_N][-1] = -1.0
    
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

    def plot_lattice(self,layer = 0, aspect = 6, show_staple_nodes = False, title: str = "Crossover Positions"):
        lattice_2D = self.lattice[:,:, layer]
        fig, ax = plt.subplots(figsize=(15,15))
        # Set Colours
        vals = np.linspace(0,1,256)
        # np.random.shuffle(vals)
        # cmap = plt.cm.colors.ListedColormap(plt.cm.gnuplot(vals))
        cmap = plt.cm.colors.ListedColormap(plt.cm.viridis(vals))
        cmap.set_bad(color='black')
        # Create plot
        im = ax.imshow(lattice_2D, cmap=cmap,aspect = 2)
          
        # Add text overlay
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
                    ax.text(j, i, text,ha="center", va="center", color="w", fontsize=20)
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
                    ax.text(j, i, text,ha="center", va="center", color="w", fontsize=20)

        # Make it look a little prettier
        ax = plt.gca()
        ax.set_aspect(aspect)
        ax.invert_yaxis()
        ax.set_title(title)
        ax.set_xlabel("No. of nucleotides")
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.set_ylabel("Scaffold row strands")


        # Add For Staples Colorbar
        # def discrete_matshow(data, cmap, name, shrink_val: float, max_data, min_data, space_between_ints):
        #     #get discrete colormap

        #     cmap = plt.get_cmap('viridis', int(max_data-min_data+1))        
        #     cmap.set_bad(color='black')
        #         # set limits .5 outside true range
        #     mat = plt.matshow(data,cmap=cmap,vmin = min_data-.5, vmax = max_data+.5)
        #     #tell the colorbar to tick at integers
        #     cax = ax.figure.colorbar(
        #         mat, 
        #         fraction=0.046, 
        #         pad=0.04, 
        #         shrink = shrink_val, 
        #         ticks=np.arange(min_data,max_data+1,space_between_ints))
            
        #     cax.ax.set_ylabel(name, rotation = "90")

        # if layer == 0: # NOTE TO SELF doesn't work for non straight edges shapes...
        #     max_data = np.max(lattice_2D)
        #     min_data = np.min(lattice_2D)
        #     discrete_matshow(lattice_2D, cmap, "Staple ID", 0.43, max_data, min_data, 1)
        # elif layer == 3:
        #     max_data = 1930.0 #np.max(data)
        #     min_data = 0.0 #np.min(data)
        #     discrete_matshow(lattice_2D, cmap, "Scaffold Nucleotide ID", 0.95, max_data, min_data, max_data)
            
        # Show
        plt.rcParams['font.size'] = '22'  
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
        # assert len(staples) <= 2
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

    def generate_origami(self):
        """ Converts `self.lattice` to staples using `self.scaffold_obj` returning `OrigamiContainer` object

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
        
        return OrigamiContainer(staple_strands = staples, staple_base_class = self)
    
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
        self.plot_lattice(layer = 0, show_staple_nodes=False, title="Initialised Staple Framework waith crossovers")
        self.domain_size = domain_size
        print(f"StapleAlgorithm1: Generating side staples {self}")
        self.generate_side_staples()
        print(f"StapleAlgorithm1: Generating inside staples of single domain only {self}")
        # self.generate_inside_staples()
    
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

        self.plot_lattice(layer = 0, show_staple_nodes=False, title="Initialised Staple Framework with crossovers")

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
        # Set up
        self.i_left = 0
        if self.info_dataframe['start side'][0] == "left":
            self.i_right = self.info_dataframe['5p'][0]
        else:
            logger.kill(f"Scaffold route travelling in this direction is currently unsupported")
        
        # Calculate
        self.i_q1 = int( round ( 1/4 * (self.i_right-self.i_left) , 0 ) )
        self.i_q2 = int( round ( 2/4 * (self.i_right-self.i_left) , 0 ) )
        self.i_q3 = int( round ( 3/4 * (self.i_right-self.i_left) , 0 ) )

    def find_closest_crossover_indices(self):
        """ Ensures 5 potential crossover points all occur where  """

        def _find_crossover(guess_idx: int, row: int) -> int:
            """ Oscillates around guess_idx to find a nt where the
            a1 vector of the scaffold is pointing down (== -1)"""   
            modifier_list = [-1,2,-3,4,-5,6,-7,8,-9]
            i = 0
            while not self.lattice[row, guess_idx, 1] == -1:
                # modify index by oscillating around the current index +1, -2, +3... etc
                guess_idx += modifier_list[i]
                i += 1

            if not self.lattice[row, guess_idx, 1] == -1:
                logger.kill(f"Crossover point was not found near {guess_idx} on row {row}")
            else:
                logger.info(f"Crossover was found at an index of {guess_idx} on row {row}")
                return guess_idx
        
        self.i_left =   _find_crossover(self.i_left, 0)
        self.i_q1   =   _find_crossover(self.i_q1,   0)
        self.i_q3   =   _find_crossover(self.i_q3,   1)
        self.i_q2   =   int(round((self.i_q1+self.i_q3)/2,0))
        # self.i_q2   =   _find_crossover(self.i_q2,   1) no crossovers occur at this point
        self.i_right =  _find_crossover(self.i_right,1)

        def duplicate_check(list_of_integers: list) -> bool:
            """ Check for duplicates"""
            for item in list_of_integers:
                if list_of_integers.count(item) > 1:
                    return True
                return False

        if duplicate_check([self.i_left, self.i_q1,self.i_q2,self.i_q3,self.i_right]):
            logger.kill(f"A unique set of crossover points was not found")
        else:
            logger.info(f"A set of suitable crossovers has been found")

    
    def generate_staples(self):
        """ Generates staples in two cycles: 1) left half, 2) right half """
        logger.info("StapleAlgorithm2: Generating staples")
        staples = []
        if self.info_dataframe['start side'][0] == "left":
            ## ADD STAPLES TO LEFT HALF
            for row in range(0,self.n_rows,2):
                if not row == self.n_rows - 1: # not the last row
                    logger.info("Staples left on row {} and {}".format(row, row+1))  
                    staple_left_1 = [
                        [self.i_q1-1, row,   'S'],
                        [self.i_left, row,   'C'],
                        [self.i_left, row+1, 'C'],
                        [self.i_q1-1, row+1, 'T'],
                    ]
                    staple_left_2 = [
                        [self.i_q2-1, row,   'S'],
                        [self.i_q1,   row,   'C'],
                        [self.i_q1,   row+1, 'C'],
                        [self.i_q2-1, row+1, 'T'],
                    ]
                else:
                    logger.info("Staples left on last row {}".format(row))  
                    staple_left_1 = [[self.i_q1-1, row, 'S'], [self.i_left, row, 'T']]
                    staple_left_2 = [[self.i_q2-1, row, 'S'], [self.i_q1, row, 'T']]
                
                staples.append(staple_left_1)
                staples.append(staple_left_2)
            
            ## ADD STAPLES TO RIGHT HALF
            for row in range(1,self.n_rows,2):
                if not row == self.n_rows - 1: # not the last row
                    logger.info("Staples right on row {} and {}".format(row, row+1))  
                    staple_right_1 = [
                        [self.i_q2, row,   'S'],
                        [self.i_q3, row,   'C'],
                        [self.i_q3, row+1, 'C'],
                        [self.i_q2, row+1, 'T'],
                    ]
                    staple_right_2 = [
                        [self.i_q3+1,  row,   'S'],
                        [self.i_right, row,   'C'],
                        [self.i_right, row+1, 'C'],
                        [self.i_q3+1,  row+1, 'T'],
                    ]
                else:
                    logger.info("Staples right on last row {}".format(row))  
                    staple_right_1 = [[self.i_q2,   row, 'S'], [self.i_q3,    row, 'T']]
                    staple_right_2 = [[self.i_q3+1, row, 'S'], [self.i_right, row, 'T']]
                
                staples.append(staple_right_1)
                staples.append(staple_right_2)
            
            # Single domain for row 0
            staple_right_1 = [[self.i_right, 0, 'S'], [self.i_q3+1, 0, 'T']]
            staple_right_2 = [[self.i_q3,    0, 'S'], [self.i_q2,   0, 'T']]
            
            staples.append(staple_right_1)
            staples.append(staple_right_2)

            self.write_staples_to_array(staples)
        else:
            logger.kill("Stapling for this direction has not yet been implemented")

class StaplingAlgorithm3(StapleBaseClass):
    """
    UNFINISHED algorithm
    This algorithm is specifically for rectangular shapes. i.e. the edges must be straight

    It goes about putting x+2 number of staples on each scaffold row. Every row will always have a
    staple crossover holding two scaffold crossovers together, hence (+2).

    Split up as x/2 from the left, x/2 from the right (shifted up one row).

    Steps:
    1. Initialise size of crossover_points list: self.cp = np.array(x+2)
    2. Calculate approximate index of staple crossover positions for 2 + x (start/end + middle) crossover points
       for staples added to the scaffold
        - L  (First nt) 
        - M  (nt at approx the 1,,,x/(x+1) mark)
        - R  (Last nucleotide)

    3. Do an oscillatory search around those points for the coordinates 
        corresponding to a crossover in the scaffold
    3. 
    """
    def __init__(self, scaffold: LatticeRoute, middle_staples: int, crossover_threshold=0.956):
        super().__init__(scaffold, crossover_threshold)
        self.middle_staples = middle_staples
        self.cp = np.zeros(self.middle_staples+2) # 2 = 1 left + 1 right
        print(self.info_dataframe)
        self.verify_route_is_applicable()
        self.find_approx_crossover_indices()
        # self.find_closest_crossover_indices()
        # self.generate_staples()


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
        # Set up
        self.i_left = 0
        if self.info_dataframe['start side'][0] == "left":
            self.i_right = self.info_dataframe['5p'][0]
        else:
            logger.kill(f"Scaffold route travelling in this direction is currently unsupported")
        
        for i in range(len(self.cp)):
            self.cp[i] = int( round ( i/(len(self.cp)-1) * (self.i_right-self.i_left), 0 ) )

        assert self.cp[0] == self.i_left, f"{self.cp[0]},{self.i_left}"
        # logger.info(f"The crossover points are {self.cp}")
        assert self.cp[-1] == self.i_right, f"{self.cp[-1]},{self.i_right}"

    def find_closest_crossover_indices(self):
        """ Ensures 5 potential crossover points all occur where  """

        def _find_crossover(guess_idx: int, row: int) -> int:
            """ Oscillates around guess_idx to find a nt where the
            a1 vector of the scaffold is pointing down (== -1)"""   
            modifier_list = [1,-2,3,-4,5,-6, 7, -8, 9, -10, 11, -12] 
            i = 0
            while not self.lattice[row, guess_idx, 1] == -1:
                # modify index by oscillating around the current index +1, -2, +3... etc
                guess_idx += modifier_list[i]
                i += 1

            if not self.lattice[row, guess_idx, 1] == -1:
                logger.kill(f"Crossover point was not found near {guess_idx} on row {row}")
            else:
                logger.info(f"Crossover was found at an index of {guess_idx} on row {row}")
                return guess_idx
        
        def duplicate_check(list_of_integers: list) -> bool:
            """ Check for duplicates"""
            for item in list_of_integers:
                if list_of_integers.count(item) > 1:
                    return True
                return False

        self.i_left =   _find_crossover(self.i_left, 0)
        self.i_q1   =   _find_crossover(self.i_q1,   0)
        self.i_q3   =   _find_crossover(self.i_q3,   1)
        # self.i_q2   =   int(round((self.i_q1+self.i_q3)/2,0))
        # self.i_q2   =   _find_crossover(self.i_q2,   1) no crossovers occur at this point
        self.i_right =  _find_crossover(self.i_right,1)

        if duplicate_check([self.i_left, self.i_q1,self.i_q2,self.i_q3,self.i_right]):
            logger.kill(f"A unique set of crossover points was not found")
        else:
            logger.info(f"A set of suitable crossovers has been found")

    
    def generate_staples(self):
        """ Generates staples in two cycles: 1) left half, 2) right half """
        logger.info("StapleAlgorithm2: Generating staples")
        staples = []
        if self.info_dataframe['start side'][0] == "left":
            ## ADD STAPLES TO LEFT HALF
            for row in range(0,self.n_rows,2):
                if not row == self.n_rows - 1: # not the last row
                    logger.info("Staples left on row {} and {}".format(row, row+1))  
                    staple_left_1 = [
                        [self.i_q1-1, row,   'S'],
                        [self.i_left, row,   'C'],
                        [self.i_left, row+1, 'C'],
                        [self.i_q1-1, row+1, 'T'],
                    ]
                    staple_left_2 = [
                        [self.i_q2-1, row,   'S'],
                        [self.i_q1,   row,   'C'],
                        [self.i_q1,   row+1, 'C'],
                        [self.i_q2-1, row+1, 'T'],
                    ]
                else:
                    logger.info("Staples left on last row {}".format(row))  
                    staple_left_1 = [[self.i_q1-1, row, 'S'], [self.i_left, row, 'T']]
                    staple_left_2 = [[self.i_q2-1, row, 'S'], [self.i_q1, row, 'T']]
                
                staples.append(staple_left_1)
                staples.append(staple_left_2)
            
            ## ADD STAPLES TO RIGHT HALF
            for row in range(1,self.n_rows,2):
                if not row == self.n_rows - 1: # not the last row
                    logger.info("Staples right on row {} and {}".format(row, row+1))  
                    staple_right_1 = [
                        [self.i_q2, row,   'S'],
                        [self.i_q3, row,   'C'],
                        [self.i_q3, row+1, 'C'],
                        [self.i_q2, row+1, 'T'],
                    ]
                    staple_right_2 = [
                        [self.i_q3+1,  row,   'S'],
                        [self.i_right, row,   'C'],
                        [self.i_right, row+1, 'C'],
                        [self.i_q3+1,  row+1, 'T'],
                    ]
                else:
                    logger.info("Staples right on last row {}".format(row))  
                    staple_right_1 = [[self.i_q2,   row, 'S'], [self.i_q3,    row, 'T']]
                    staple_right_2 = [[self.i_q3+1, row, 'S'], [self.i_right, row, 'T']]
                
                staples.append(staple_right_1)
                staples.append(staple_right_2)
            
            # Single domain for row 0
            staple_right_1 = [[self.i_right, 0, 'S'], [self.i_q3+1, 0, 'T']]
            staple_right_2 = [[self.i_q3,    0, 'S'], [self.i_q2,   0, 'T']]
            
            staples.append(staple_right_1)
            staples.append(staple_right_2)

            self.write_staples_to_array(staples)
        else:
            logger.kill("Stapling for this direction has not yet been implemented")


class OrigamiContainer:
    """ Stores staple and scaffold strand objects with the ability to modify them.

    T0D0: add methods to modify staples -> could be added through subclassing
    """
    def __init__(self, 
                staple_strands: List[Strand] = [], 
                staple_base_class: StapleBaseClass = None):

        self._staples = staple_strands
        self._baseclass = staple_base_class
        self._scaffold = self._baseclass.scaffold_obj
        self._scaffold_rows = self._baseclass.scaffold_rows
    
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
        """ Generates a plot for the scaffold and original set of staples generated in the StapleAlgorithm """
        self._baseclass.plot_lattice()

class Configuration:
    """ Class object which stores multiple strands of DNA and can export/write to file
    
    Simplified version of OrigamiContainer
    """
    def __init__(self, staples: list, scaffold: Strand):
        self.staples = staples
        self.scaffold = scaffold

    @property
    def dna_strands(self) -> List[Strand]:
        """ Returns scaffold and corresponding staple strands """
        return self.staples + [self.scaffold]

    @property
    def n_staples(self) -> int:
        return len(self.dna_strands) - 1

    def system(self, **kwargs) -> System:
        """ Returns oxDNA.System object with the scaffold and all staple strands """
        _system = System(kwargs.get('box', np.array([50., 50., 50.])))
        _system.add_strands(self.dna_strands)
        return _system
        
    def output_files(self, name: str, root_dir: str = ".", **kwargs):
        """ Generates oxDNA files (`.top` and `.conf`) for the strand configuration """
        logger.info(f"Writing oxDNA files: oxdna.{name}.top/conf ")
        return self.system(**kwargs).write_oxDNA(prefix = name, root = root_dir)


class ConfGenSplitDoubleDomains(OrigamiContainer):
    """
    Contains functions to split a double domained staple into two single
    -domained staples. Stores these in `Configuration` objects, which can be outputted
    as oxDNA and LAMMPS objects.
    """
    def __init__(self, 
                staple_strands: List[Strand] = [], 
                staple_base_class: StapleBaseClass = None):
        super().__init__(staple_strands, staple_base_class)
        logger.info("Configuration Generator initialised.")
        ## Initialise
        self.configurations = []
        self.configurations.append(Configuration(self.staples, self.scaffold))
        # Min/max x-coord of scaffold
        self.min_x, self.max_x = self.find_minmax_x_of_scaffold()

        ## Generate Configurations
        logger.info("Generating multiple configurations")
        self.generate_all_configs_1()


    def find_minmax_x_of_scaffold(self):
        """ Find min and max x coord of nt backbone position on scaffold.
        
        First collect x positions of all nt on scaffold
        """
        scaf = self._scaffold_rows[0]
        backbone_x = [nt.pos_back[0] for nt in scaf.nucleotides]
        _min = min(backbone_x)
        _max = max(backbone_x)
        return _min, _max

    def inside_staple(self, staple: Strand) -> bool:
        """ 
        Returns true if a given staple does not have any nts at the edges of the scaffold.
        
        Checks whether a staple contains a max and min x-coord of the backbone (nucleotide.pos_back[0])
        and if it does, then it means that the staple is NOT an inside staple.
        """
        backbone_x = [nt.pos_back[0] for nt in staple.nucleotides]
        
        return not any([any(x <= self.min_x for x in backbone_x), any(x >= self.max_x for x in backbone_x)])

    @staticmethod
    def double_domained(staple: Strand) -> bool:
        """ Checks whether a staple is double domained """
        # Ensure staple is double domained
        nt_dir = [nt._a3[0] for nt in staple.nucleotides] # list of -1.0 (and if double domained 1.0)
        unique_directions = set(nt_dir)
        if len(unique_directions) == 1:
            print("Staple is single domained -> skipping")
            return False
        else:
            return True

        
    @staticmethod
    def split_staple(staple_to_split: Strand):
        """ Splits a double domained staple into two single domained staples.

        Finds the index of the nt where the backbone site changes by 1 unit
        in the y direction (i think). Then copies the nt up to that point and 
        again copies the nt after that point, forming two new oxDNA strand objects.

        Returns (oxDNA.Strand, oxDNA.Strand):
            a tuple of two single domained staples 
        """
        # We want to ensure the direction has flipped... so we find the a3 direction
        # This will be either 1.0 or -1.0
        # If it starts at 1.0, we are looking for the nucleotide where the a3 direction is -1.0
        # and vice versa

        # Find index of first nt in the second domain of the staple
        i = 0
        hello = 0
        target = -staple_to_split.nucleotides[i]._a3[0] # negative of result, i.e. 1 -> -1, -1 -> 1

        while not hello == target:
            i += 1
            hello = staple_to_split.nucleotides[i]._a3[0]
        
        # if i == 0:
        #     logger.info(f"{[(nt._a1[0], nt._a2[0], nt._a3[0]) for nt in staple_to_split.nucleotides]}")
        
        logger.info(f"THIS STAPLE IS BEING CUT AT INDEX {i} and its length is {len(staple_to_split.nucleotides)}")

        # Define new single domain strands
        strand_1 = Strand(nucleotides=staple_to_split.nucleotides[:i])
        strand_2 = Strand(nucleotides=staple_to_split.nucleotides[i:])

        return strand_1, strand_2

    def generate_all_configs_1(self):
        """ Generates configurations with one inner staple converted to two single domains

        Yields additions to list of configuration objects (which can be written to file)
        """
        for i, staple in enumerate(self.staples):
            if self.inside_staple(staple) and self.double_domained(staple):

                # Copy list of staples and remove staple we want to split
                new_staple_set = deepcopy(self.staples)
                logger.info(f"Replacing staple {i} with two single domain staples")
                
                # Split double domain staple to two single domain staples
                new_staple_1, new_staple_2 = self.split_staple(staple)

                # Add new staples to list             
                new_staple_set.append(new_staple_1)
                new_staple_set.append(new_staple_2)
                del new_staple_set[i]

                # Create configuration and append to list of configurations
                new_configuration = Configuration(new_staple_set, self.scaffold)
                assert new_configuration.n_staples == len(new_staple_set)
                assert len(new_staple_set) == len(self.staples) + 1
                print("Generating a new configuration of staples")
                self.configurations.append(new_configuration)

    def write_to_file(self, name: "str", root = "."):
        tot = len(self.configurations)
        for i, conf in enumerate(self.configurations):
            name_str = name + "-" + str(i)
            logger.info(f"Writing files {i+1}/{tot}")
            conf.output_files(name_str, root)
    

class ConfGenAddUnpairedNuc(OrigamiContainer):
    """ Configuration Generator: Add Unpaired Nucleotides to Scaffold/Staple crossovers

    Contains functions to find crossover of DNA strand
    and add x no of nucleotides to that staple

    Parameters:
        staple_strands - staple strands forming the origami, as a list of oxDNA strand objects
        staple_base_class - the object created through all `StaplingAlgorithm`s
    """
    def __init__(self,
                staple_strands: List[Strand] = [], 
                staple_base_class: StapleBaseClass = None):
        super().__init__(staple_strands, staple_base_class)
        logger.info(f"Configuration Generator: 'Add Unpaired Nucleotides' initialised.")
        self.configuration = Configuration(self.staples, self.scaffold)
        self._new_configuration = None

    @property
    def new_configuration(self) -> Configuration:
        if self._new_configuration is not None:
            return self._new_configuration
        else:
            logger.error(f"You must run `.add_to_strand/scaffold` first")

    @staticmethod
    def add_nt_to_crossovers(strand: Strand) -> Strand:
        """ Adds (currently a single) nt to the given oxDNA strand object at all crossover point """
        direction = [nt._a3[0] for nt in strand.nucleotides]

        # Find all indices of crossovers in strand
        crossover_index = np.array([index for index, value in enumerate(direction[:-2]) if value != direction[index+1]])
        for index in crossover_index:
            assert direction[index+1] != direction[index]
        crossover_index_pairs = np.stack((crossover_index, crossover_index + 1), axis = 1)

        if len(crossover_index_pairs) == 0:
            logger.info("Not adding nt(s) to staple as it is single domained (no crossover exists)")   

        # Add (1) nucleotide to all crossovers
        else:
            for added_nt, pair_indicies in enumerate(crossover_index_pairs):
                # update list/index of nucleotides for each loop - accounts for inserted nt in previous loops
                nucleotides = strand.nucleotides
                [n1_idx, n2_idx] = pair_indicies + added_nt

                # position COM of new base
                new_pos_com = (nucleotides[n1_idx].pos_com + nucleotides[n2_idx].pos_com)/2
                new_pos_com += nucleotides[n1_idx]._a3 * 0.5  # shift horizontally
                # new_pos_com -= nucleotides[n1_idx]._a2 * POS_BASE * 0.3 # shift COM such that backbone is inbetween n1 and n2
                new_pos_com -= nucleotides[n1_idx]._a1 * 0.4

                # define orientation/tilting vectors
                new_a1 = -nucleotides[n1_idx]._a1
                new_a3 = nucleotides[n1_idx]._a3
                
                # create new nucleotide
                new_nt = Nucleotide(
                    random.choice(['A', 'T', 'C', 'G']),
                    new_pos_com,
                    new_a1,
                    new_a3
                    )
                # insert nucleotide before given index
                strand.nucleotides.insert(n2_idx, new_nt) 

        return strand

    @staticmethod
    def add_2nt_to_crossovers(strand: Strand) -> Strand:
        """ Adds (currently a single) nt to the given oxDNA strand object at all crossover point """
        direction = [nt._a3[0] for nt in strand.nucleotides]

        # Find all indices of crossovers in strand
        crossover_index = np.array([index for index, value in enumerate(direction[:-2]) if value != direction[index+1]])
        for index in crossover_index:
            assert direction[index+1] != direction[index]
        crossover_index_pairs = np.stack((crossover_index, crossover_index + 1), axis = 1)

        # check for single domain staples
        if len(crossover_index_pairs) == 0:
            logger.info("Not adding nt(s) to staple as it is single domained (no crossover exists)")   

        # Add 2 nucleotides to all crossovers
        else:
            for added_nt, pair_indicies in enumerate(crossover_index_pairs):
                # update list/index of nucleotides for each loop - accounts for inserted nt in previous loops
                nucleotides = strand.nucleotides
                [n1_idx, n2_idx] = pair_indicies + added_nt*2

                # shift both bases in the same a3 direction (n1._a3)
                new_1_pos_com = nucleotides[n1_idx].pos_com + 0.5 * nucleotides[n1_idx]._a3
                new_2_pos_com = nucleotides[n2_idx].pos_com + 0.5 * nucleotides[n1_idx]._a3
                
                # shift COM towards center of crossover (y) direction
                new_1_pos_com -= nucleotides[n1_idx]._a1 * 0.1
                new_2_pos_com -= nucleotides[n2_idx]._a1 * 0.1
                
                # orientation and tilting of backbone
                new_1_a1 = nucleotides[n1_idx]._a1
                new_1_a3 = nucleotides[n1_idx]._a3 
                new_2_a1 = nucleotides[n2_idx]._a1
                new_2_a3 = nucleotides[n2_idx]._a3 
                
                # create new nucleotides 
                new_1_nt = Nucleotide(
                    random.choice(['A', 'T', 'C', 'G']),
                    new_1_pos_com,
                    new_1_a1,
                    new_1_a3
                    )
                new_2_nt = Nucleotide(
                    random.choice(['A', 'T', 'C', 'G']),
                    new_2_pos_com,
                    new_2_a1,
                    new_2_a3
                    )
                # insert nucleotides before given index 
                strand.nucleotides.insert(n2_idx, new_2_nt)
                strand.nucleotides.insert(n2_idx, new_1_nt) 

        return strand

    def initialise_new_conf(self):
        """ Sets old configuration = new configuration """
        self._new_configuration = self.configuration

    def add_to_scaffold(self, add_x_nt: int):
        """ Adds *x* nucleotide(s) to every scaffold crossover 
        
        This function works by retrieving current scaffold from self.new_configuration, 
        and adds one nucleotide to all crossovers. 

        Parameter:
            add_x_nt - number of nucleotides to add to scaffold crossover

        Yields:
            updated self._new_configuration

        """
        if self._new_configuration is None:
            self.initialise_new_conf()
        
        if add_x_nt == 1:
            logger.info("Adding *1* nucleotide to every scaffold crossover.")
            new_scaffold = self.add_nt_to_crossovers(self._new_configuration.scaffold)
        elif add_x_nt == 2:
            logger.info("Adding *2* nucleotide to every scaffold crossover.")
            new_scaffold = self.add_2nt_to_crossovers(self._new_configuration.scaffold)
        else:
            logger.error("Addition of only 1 or 2 nt is currently supported!")
        
        self._new_configuration.scaffold = new_scaffold

    def add_to_all_staples(self, add_x_nt: int):
        """ Adds *x* nucleotide(s) to every staple crossover 

        This function works by retrieving current scaffold from self.new_configuration, 
        and adds one nucleotide to all crossovers. 

        Parameter:
            add_x_nt - number of nucleotides to add to scaffold crossover

        Yields:
            updated self._new_configuration

        """
        if self._new_configuration is None:
            self.initialise_new_conf()
        
        if add_x_nt == 1:
            logger.info("Adding *1* nucleotide to every staple at its crossover.")
            staples = [self.add_nt_to_crossovers(staple) for staple in self._new_configuration.staples]
        elif add_x_nt == 2:
            logger.info("Adding *2* nucleotides to every staple at its crossover.")
            staples = [self.add_2nt_to_crossovers(staple) for staple in self._new_configuration.staples]
        else:
            logger.error("Addition of only 1 or 2 nt is currently supported!")
        
        self._new_configuration.staples = staples

        return

## Copied from protocols/lattice-route/DNA_snake.py
def generate_scaffold(polygon_vertices: np.ndarray, title: str = "generate_scaffold()", bp_per_turn: float = 10.45) -> LatticeRoute:
    print('Running Generate Function')
    print(f'Creating polygon from {len(polygon_vertices)} vertices')
    polygon = BoundaryPolygon(polygon_vertices)
    print("Constructing scaffolding lattice")
    # print(f"{title}: ...constructing scaffolding lattice...")
    lattice = polygon.dna_snake(start_side="left", grid_size = [0.34, 2], bp_per_turn = bp_per_turn)
    # print(f"{title}: ...calculating route.")
    print(f"Calculating Route for lattice with {lattice.bp_per_turn} helical twist")
    route = lattice.route()
    print("Plotting Route")
    route.plot()
    return route

## Protocol functions
def staple_1_and_write_to_file(route: LatticeRoute, name_of_file: str, domain_size = 15):
    logger.info("Class StaplingAlgorithm1: Generating side staples")
    with_staples_1 = StaplingAlgorithm1(route, domain_size = domain_size)

    with_staples_1.plot_lattice(title=name_of_file)

    logger.info("staple_1_and_write_to_file(): Adding staples to a container...")
    container_1 = with_staples_1.generate_origami()
    
    logger.info("staple_1_and_write_to_file(): Adding staples to an oxDNA system...")
    system_1 = container_1.system()
    
    logger.info("staple_1_and_write_to_file(): Writing `.top` and `.conf` files")
    system_1.write_oxDNA(prefix = name_of_file)
    return with_staples_1, container_1

def staple_2_and_write_to_file(route: LatticeRoute, name_of_file: str):
    logger.info("Class StaplingAlgorithm2: Generating 4 Columns of Staples")
    staple_2 = StaplingAlgorithm2(route)
    
    staple_2.plot_lattice(title=name_of_file)
    
    logger.info("staple_1_and_write_to_file(): Adding staples to a container...")
    container_2 = staple_2.generate_origami()

    logger.info("staple_1_and_write_to_file(): Adding staples to an oxDNA system...")
    # system_2 = container_2.system()

    logger.info("staple_1_and_write_to_file(): Writing `.top` and `.conf` files")
    # system_2.write_oxDNA(name)
    return staple_2, container_2

def plot_staples(staples: OrigamiContainer):
    if type(staples) == OrigamiContainer:
        staples = staples.base_class

    print("Plotting Staples with start, crossover and terminal points labelled")
    staples.plot_lattice(layer=0, show_staple_nodes=True, title="Layer 0 with start, crossover and terminal points labelled")
    print("Plotting Staples with scaffold crossover points labelled")
    staples.plot_lattice(layer=0, show_staple_nodes=False, title="Layer 0 with scaffold crossover points labelled")
    staples.plot_lattice(layer=1, show_staple_nodes=False, title="Layer 1 with scaffold crossover points labelled")
    staples.plot_lattice(layer=2, show_staple_nodes=False, title="Layer 1 with scaffold crossover points labelled")
    staples.plot_lattice(layer=3, show_staple_nodes=False, title="Layer 1 with scaffold crossover points labelled")

def param_study_0002():
    """ A parameter study looking at how the change of a double domain to single domain affects the structure """
    square = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])*np.array([3.5,2.5,1])
    rectangle = square*[8,5,1]
    route = generate_scaffold(rectangle)
    
    staple_2 = StaplingAlgorithm2(route)
    staple_2.plot_lattice()
    container_2 = staple_2.generate_origami()
    system_2 = container_2.system()
    # system_2.write_oxDNA("half")

    configgen = ConfGenSplitDoubleDomains(staple_strands = container_2.staples, staple_base_class = staple_2)
    ROOT = "/".join(path.abspath(__file__).split("/")[:-1])
    print(ROOT)
    configgen.write_to_file(name = "batch1", root = ROOT) 
    return staple_2, container_2

def param_study_0002_but_square():
    """ A parameter study looking at how the change of a double domain to single domain affects the structure """
    
    stapled_scaffold, staple_container = generate_origami_square(width = 4)

    configgen = ConfGenSplitDoubleDomains(staple_strands = staple_container.staples, staple_base_class = stapled_scaffold)
    ROOT = "/".join(path.abspath(__file__).split("/")[:-1])
    print(ROOT)
    configgen.write_to_file(name = "square-removed_crossovers", root = ROOT) 
    return stapled_scaffold, staple_container

def generate_origami_square(width: float=8, name: str="Stapled Scaffold Schematic"):
    ### Pick A Shape
    square = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])*np.array([3.5,12.5,1])

    ### Pick An Algorithm
    # route = generate_scaffold(stacked_I*[12,8,8])
    route = generate_scaffold(square*[width,1,1])
    scaffold_system = route.system()
    scaffold_system.write_oxDNA("scaffold_system")
    # staple, container = staple_1_and_write_to_file(route, "square25", domain_size=25)
    # plot_staples(container)

    staple_2, container_2 = staple_2_and_write_to_file(route, "square")
    origami_system = container_2.system()
    origami_system.write_oxDNA("origami_system")
    plot_staples(container_2)
    return staple_2, container_2

def generate_origami_stacked_I(width: float=8, name: str="Stapled Scaffold Schematic"):
    ### Pick A Shape
    hourglass = np.array([[0.,0.,0.],[4,6.,0.],[0,12,0],[12,12,0],[8,6.,0.],[12,0.,0.]])
    stacked_I = np.array([
    [0.,0.,0.],[3.,0.,0.],[3.,1.,0.],[2.,1.,0.], [2.,2.,0.],[3.,2.,0.],
    [3.,3.,0.],[2.,3.,0.],[2.,4.,0.],[3.,4.,0.],[3.,5.,0.],[0.,5.,0.],[0.,4.,0.],[1.,4.,0.],
    [1.,3.,0.],[0.,3.,0.], [0.,2.,0.],[1.,2.,0.],[1.,1.,0.],[0.,1.,0.]
    ])
    square = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])*np.array([3.5,12.5,1])
    triangle = np.array([[0,0,0],[5,10,0],[10,0,0]])
    trapREV = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]])

    ### Pick An Algorithm
    route = generate_scaffold(stacked_I*[12,8,8])
    scaffold_system = route.system()
    scaffold_system.write_oxDNA("scaffold_system_stackedI")
    staple, container = staple_1_and_write_to_file(route, "stacked_I", domain_size=25)
    plot_staples(container)

    # staple_2, container_2 = staple_2_and_write_to_file(route, "square")
    # origami_system = container_2.system()
    # origami_system.write_oxDNA("origami_system")
    # plot_staples(container_2)
    # return staple_2, container_2

def param_study_0003():
    """ A study looking at the same staple architecture with varying aspect ratio"""
    for i in np.linspace(2,10,9):
        width = round(i,2)
        shape_name = "square-width-" + str(width)
        staple_2, container_2 = generate_origami_square(width, shape_name)

def param_study_0006():
    """ Adding 1 nt to all scaffold crossovers """
    # Make Square-width-4.0
    stapled_scaffold, staple_container = generate_origami_square(width = 4)
    configuration_container = ConfGenAddUnpairedNuc(
        staple_strands = staple_container.staples,
        staple_base_class = stapled_scaffold
        )
    
    origami_0 = deepcopy(configuration_container)
    origami_1 = deepcopy(configuration_container)
    origami_2 = deepcopy(configuration_container)
    
    origami_0.configuration.output_files("scaffold_0_staples_0")
    origami_0.add_to_all_staples(1)
    origami_0.configuration.output_files("scaffold_0_staples_1")
    origami_0 = deepcopy(configuration_container)
    origami_0.add_to_all_staples(2)
    origami_0.configuration.output_files("scaffold_0_staples_2")    

    origami_1.add_to_scaffold(1)
    origami_1.configuration.output_files("scaffold_1_staples_0")
    origami_1.add_to_all_staples(1)
    origami_1.configuration.output_files("scaffold_1_staples_1")
    origami_1 = deepcopy(configuration_container)
    origami_1.add_to_scaffold(1)
    origami_1.add_to_all_staples(2)
    origami_1.configuration.output_files("scaffold_1_staples_2")    


    origami_2.add_to_scaffold(2)
    origami_2.configuration.output_files("scaffold_2_staples_0")
    origami_2.add_to_all_staples(1)
    origami_2.configuration.output_files("scaffold_2_staples_1")
    origami_2 = deepcopy(configuration_container)
    origami_2.add_to_scaffold(2)
    origami_2.add_to_all_staples(2)
    origami_2.configuration.output_files("scaffold_2_staples_2")    

    return configuration_container
    #configuration_container.new_configuration

def param_study_000X():
    """ DOESN'T WORK! Looking at the effect of twist on a small origami """
    square = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])*np.array([3.5,2.5,1])*[3,5,1]
    for helical_twist in np.arange(10.45,10.65,0.05):
        twist = round(helical_twist,2)
        # Make Scaffold
        route = generate_scaffold(square, bp_per_turn = twist)
        # Staple Scaffold
        stapled_route = StaplingAlgorithm2(route)
        # Plot
        title = f"Stapled Scaffold Schematic (helical twist: {twist}"
        stapled_route.plot_lattice(title=title)
        # Export as oxDNA
        name = f"twist-{twist:.2f}"
        container = stapled_route.generate_origami()
        system = container.system()
        system.write_oxDNA(name)

    return route

def staple_3():
    square = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])*np.array([3.5,2.5,1])
    rectangle = square*[8,5,1]
    route = generate_scaffold(rectangle)
    StaplingAlgorithm3(route, middle_staples=2)

def for_presentation():
    # Define Shape
    square = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])*np.array([3.5,12.5,1])
    polygon = BoundaryPolygon(vertices = square*[4,1,1]) # square width 4.0
    
    # Generate Scaffold
    lattice = polygon.dna_snake(grid_size = [0.34, 2])
    scaffold = lattice.route()
    scaffold_system = scaffold.system()
    scaffold_system.write_oxDNA("scaffold_only")

    # Generate Staples
    staples = StaplingAlgorithm2(scaffold)
    origami = staples.generate_origami()
    origami_system = origami.system()
    origami_system.write_oxDNA("full_origami")

    # Generate Additional Configurations
    origami_container = ConfGenAddUnpairedNuc(origami.staples, origami.base_class)
    new_origami_1 = deepcopy(origami_container)
    new_origami_2 = deepcopy(origami_container)
    new_origami_1.add_to_scaffold(1)
    new_origami_2.add_to_scaffold(2)
    new_origami_1.configuration.output_files("extra_nt_scaffold")
    new_origami_2.configuration.output_files("extra_2nt_scaffold")

if __name__ == "__main__":
    param_study_0002()

    # param_study_0006()
    # generate_origami_stacked_I()
    # param_study_0002_but_square()

    # route = generate_scaffold(stacked_I*17)
    # system, container = staple_and_write_to_file(route, "stacked_I")
    # plot_staples(container)

    # route = generate_scaffold(hourglass*10)
    # system, container = staple_and_write_to_file(route, "hourglass")
    # plot_staples(container)

    # route = generate_scaffold(square*10)
    # system, container = staple_and_write_to_file(route, "square")
    # plot_staples(container)

    # route = generate_scaffold(trapREV*17)
    # system, container = staple_and_write_to_file(route, "trapezium")
    # plot_staples(container)

    # route = generate_scaffold(square*2)

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