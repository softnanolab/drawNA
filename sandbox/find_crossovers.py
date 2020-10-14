import itertools
from drawNA.lattice.utils import find_crossover_locations
import pandas as pd
from typing import List
from drawNA.lattice import Lattice, LatticeNode, LatticeRoute

def calculate_row_size(*half_turn_sizes: int):
    """
    Function to calculate the row size given a number of integer arguments.
    Returns: (sum - 1) and (list of half turn sizes)

    Arguments:
    half_turn_sizes - integers representing no. of nucleotides in a half turn
    """
    row_size = 0
    for half_turn_size in half_turn_sizes:
        row_size += half_turn_size
    
    # print("Half Turns at:",list(half_turn_sizes))
    row_size -= 1
    return row_size, list(half_turn_sizes)

def get_row_sizes(max_row_size, columns, bp_per_turn: float = 10.45):
    """
    Returns list and dictionary of row sizes:
        list - all possible sizes of rows as a sum of row sizes for the given number of columns
        dictionary - keys formed of values in list & values are half turn sizes which sum up to value of key

    Arguments:
        max_row_size - maximum number of lattice sites in a single row
        columns - number of pairs of crossovers on each row

    Note all key's are given as Python index values, i.e. they start from 0

    """
    half_turn_sizes = find_crossover_locations(max_row_size, bp_per_turn, origin = 1, )[1:]
    crossover_combinations = itertools.combinations_with_replacement(half_turn_sizes, columns)
    row_sizes_and_turns = list(itertools.starmap(calculate_row_size, crossover_combinations))
    
    # Generate dictionary data structure
    rowsize_dict = dict()
    for item in row_sizes_and_turns:
        key = item[0]
        value = item[1]
        rowsize_dict.setdefault(key,[]).append(value)
    
    # Generate list of only row sizes + remove duplicates + order
    rowsize_list = sorted(set([item[0] for item in row_sizes_and_turns]))
    
    return rowsize_list, rowsize_dict

rowsize_list, rowsize_dict = get_row_sizes(80, columns = 3)

df = pd.DataFrame.from_dict(rowsize_dict, orient = 'index')
df.sort_index(ascending=True, inplace=True)
# df.transpose()
display(df)


class DNAScaffold(Lattice):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        """
        Subclass of Lattice. Inherits all atributes and methods from this parent class.

        When a DNASnake object is generated, it refers to a DNA Scaffold 
        with a maximum of 2 crossovers per "row" in the lattice
        """
        self.quantised_array = self.quantise_rows(self.intersected_array, self.poss_cross, self.padding)
        # self.straightened_array = self.straighten_edges(self.quantised_array, self.straightening_factor)


        # # Final Lattice and crossovers
        # self.final_array = ?????
        # self.crossover_array = self.get_crossovers(self.final_array, self.start_side, self.poss_cross)

        # # Lattice and Crossovers as coordinates
        # self.final_coords = self.array_to_coords(self.final_array)
        # self.crossover_coords = self.array_to_coords(self.crossover_array)


        """
        Evening thoughts 27/08/20
        I really think that trying to get a way to sift through the possible crossover pairs
        (when there are multiple) for a given row size, will be hard if we just do a blanket
        "find crossovers by overlaying onto the quantised lattice" and then try to figure out an
        algorithm from there....

        Instead i think we need to create a way of:
        1. find crossover sizes for row
        2. choose the first one, if there is a second, store it for later
        3. THEN make a crossover lattice (for plotting)
        4. Look at the lattice and figure out if we need any other rules, e.g. doubling up rows or smthng
        5. Then cycle through rows and add (or even generate) vertices -> `LatticeNodes` for that row.
            Where the first set will be added to one list, the second set to another list (in reverse? idk)
            Then make a mega `vertex` list
        6. Probs gotta think about the double crossover (i.e. at the top/bottom when we want to change direction entirely). That row is gonna have to be a "whole turn" crossover, will be interesting to implement this.
        7. Once we've figured that out, we come back to looking at multiple ways of routing the scaffold 
        through the lattice
        """

    def route(self, crossovers = None, *args, **kwargs) -> List[LatticeNode]:
        """
        Generate DNASnake scaffold route, returns list of LatticeNode objects
        
        Arguments:
            crossovers - binary array or list of coordinates of crossovers
            *args & **kwargs - correspond to those of the `LatticeRoute` class
        
        Note:
            Different subclasses will use `Lattice` in a different way.
            Moreover, they will probably have a different routing algorithm
        """
        if crossovers is None:
            coords = self.crossover_coords
        elif np.shape(crossovers)[1] not in [2,3]: # if array
            coords = self.array_to_coords(crossovers)
        else:
            coords = crossovers

        # add 3rd column (z = 0)
        shape = np.shape(coords)
        if shape[1] != 3:
            coords = np.c_[coords, np.zeros(shape[0])]

        crossovers_per_row = 2
        lattice_rows = int(coords[:, 1].max() + 1)  # coords counts from 0
        vertices_in_route = int(crossovers_per_row * lattice_rows)
        vertex_list = np.zeros((vertices_in_route, 3))

        # Find final crossover from left or right and make it a node
        for row in range(0, lattice_rows):
            vertex_index_L = bisect_left(coords[:, 1], row)
            vertex_index_R = bisect_right(coords[:, 1], row) - 1
            if row % 2 == 0:  # if even
                vertex_list[row * 2] = coords[vertex_index_L]
                vertex_list[row * 2 + 1] = coords[vertex_index_R]
            else:  # if odd
                vertex_list[row * 2] = coords[vertex_index_R]
                vertex_list[row * 2 + 1] = coords[vertex_index_L]

        # print(vertex_list)

        node_list = [LatticeNode(i) for i in vertex_list]
        return LatticeRoute(node_list, *args, **kwargs)