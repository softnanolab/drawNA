from typing import List

import numpy as np
from shapely.geometry import MultiPoint
from shapely import geometry

import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

# from operator import itemgetter
from copy import deepcopy
import itertools

from .node import LatticeNode
from .route import LatticeRoute

from ..oxdna.strand import POS_STACK, FENE_LENGTH
from .utils import find_closest, find_crossover_locations
from ..oxdna.utils import round_to_multiple

from bisect import bisect_left, bisect_right

def modify_lattice_row(grid: np.ndarray,
                       difference: np.ndarray or List,
                       change_side: str = "yes"):
  
    """
    Modifies the lattice sites for a given row in the `grid` by 
    adding or removing the number of sites equal to the value which 
    is stored in `difference` for that row

    Arguments:
        grid --- (n, x) numpy array of 1's and 0's
        difference --- (n, 1) numpy array
        change_side --- ("Yes") otherwise "onlyleft" or "onlyright"
        arguments can be given, where adjustments will only be made
        on the prescribed side
    
    e.g. a difference of -5 or +11, remove (1 -> 0) or add (0 -> 1) 
    to the lattice row, respectively. Each removal/addition occurs on 
    the opposite side of the lattice than the previous removal/addition
    

    """
    # Alternate between left and right
    sides = {
        "left": "right",
        "right": "left",
        "onlyleft": "left",
        "onlyright": "right"
    }
    change_side = change_side.lower()
    message = "change_side can only have the following values: \n'yes', \n'onlyleft', \n'onlyright'"
    assert change_side in ("yes", "onlyleft", "onlyright"), message

    if not isinstance(difference, np.ndarray):
        try:
            difference = np.array([difference])
        except TypeError:
            raise TypeError(f"Difference needs to be given as a list, not {type(difference)}")
    
    if len(np.shape(grid)) != 2:
        try:
            grid = np.array([grid])
            assert len(np.shape(grid)) == 2
        except TypeError:
            raise TypeError("The value hasn't been given as a list or a np.ndarray")

    # For the difference required to reach a certain no. of nucleotides on each row
    for row, value in enumerate(difference):
        if change_side == "yes":
            side = "left"  # begin at the left side
        else:
            side = sides[change_side]

        while value != 0:
            # find location first/last 1 in the grid
            # i.e. one_right = x-coord (x location of the last 1 from the right)
            one_right = int(np.argwhere(grid[row])[:, 0].max())
            one_left = int(np.argwhere(grid[row])[:, 0].min())

            if np.sum(grid[row]) == 0:
                value = 0

            if value > 0:  #  add 1's to edges
                if side == "left":
                    grid[row, one_left - 1] = 1
                else:
                    grid[row, one_right + 1] = 1
                value -= 1
            else:  # value < 0 i.e. remove 1's from edges
                if side == "left":
                    grid[row, one_left] = 0
                else:
                    grid[row, one_right] = 0
                value += 1

            if change_side == "yes":
                side = sides[side]  # change left -> right OR right -> left

    return grid


def side(row_no: int, R: bool = None):
    """
    Returns either 'onlyleft' or 'onlyright' depending on 
    which row we are modifying, i.e. row0 will always be modified
    at its end so if R = True, then right is the end of row0,
    row1 will always be modified at the opposite side to row0
    """
    dict = {
        (0, None): "onlyleft",
        (0, True): "onlyright",
        (1, None): "onlyright",
        (1, True): "onlyleft"
    }
    return dict[(row_no, R)]


def all_same(items) -> bool:
    """Checks if all elements of a list are equal"""
    return all(x == items[0] for x in items)


class Lattice:
    def __init__(self,
                 polygon_vertices: np.ndarray,
                 grid_size: List[float] = [POS_STACK, 1.00],
                 bp_per_turn: float = 10.45,
                 straightening_factor: int = 5,
                 start_side="left"):
        """
        Lattice class forms a set of points given a polygon where a DNA Origami Scaffold
        can be laid.

        Arguments:
            polygon_vertices - a set of vertices given in order, which form a closed shape
            grid_size - the [x, y] distance between each site in the coordinate grid
                (default = [0.34, 1.00])
            bp_per_turn - the number of basepairs per 360 turn
                (default: 10.45) 
            start_side - which side the DNA scaffold will start from "left" or "right"

        Note:
            Both `lattice` and `grid` represent sites where scaffold can be laid
                `lattice`   a coordinate system
                `grid`      an array of 1's and 0's
            
        """
        self.polygon_array = polygon_vertices.astype(dtype=np.float64)
        self.polygon = geometry.Polygon(polygon_vertices)
        self.grid_size = grid_size
        self.x_spacing = grid_size[0]
        self.y_spacing = grid_size[1]
        self.bp_per_turn = bp_per_turn
        self.straightening_factor = straightening_factor
        self.padding = 40
        self.start_side = start_side.lower()

        assert self.start_side in ["left","right"], "start_side must be: 'left' or 'right'"

        self.intersected_coords = self.intersect_polygon()
        self.intersected_array = self.coords_to_array(self.intersected_coords)
        
        self.crossover_coords = None

    @property
    def poss_cross(self):
        """Find possible crossover locations"""
        max_width = np.shape(self.intersected_array)[1]
        return find_crossover_locations(max_width + self.padding)

    def intersect_polygon(self):  # Returns coordinates
        """
        Creates grid with dimensions contained in self.gridsize;
        Creates an intersection with the polygon;
        
        Returned value: set of coordinates detailing only the grid
        points which were overlaid on the polygon
        """
        # Calculate bounds of the grid around the polygon
        x_min, y_min, x_max, y_max = self.polygon.bounds
        x_min = round_to_multiple(x_min, self.x_spacing)
        x_max = round_to_multiple(x_max, self.x_spacing)
        y_min = round_to_multiple(y_min, self.y_spacing)
        y_max = round_to_multiple(y_max, self.y_spacing)

        # Initialise large grid of points spaced using x & y spacing
        # This is to "quantise" our grid
        x_grid = np.linspace(x_min, x_max,
                             int((x_max - x_min) / self.x_spacing) + 2)
        y_grid = np.linspace(y_min, y_max,
                             int((y_max - y_min) / self.y_spacing) + 1)
        bounds_grid = np.transpose(
            [np.tile(x_grid, len(y_grid)),
             np.repeat(y_grid, len(x_grid))])
        points = MultiPoint(bounds_grid)

        # Shapely intersect grid with polygon and store lattice coordinates in "lattice"
        # (Shapely isn't well documented enough and there might have been a better way
        # to do this)
        result = points.intersection(self.polygon)
        lattice = np.array(result.__geo_interface__["coordinates"])

        # Normalise values to integer values, e.g. 0.34 -> 1
        lattice[:, 0] /= self.x_spacing
        lattice[:, 1] /= self.y_spacing
        lattice = np.around(lattice).astype(int)

        # Move grid to begin at [0,0,(0)]
        x_min = lattice[:, 0].min()
        y_min = lattice[:, 1].min()
        lattice[:, 0] -= x_min
        lattice[:, 1] -= y_min
        return lattice

    @staticmethod
    def coords_to_array(coords):  # Returns array
        """ 
        Make binary array representing the lattice sites
        --> 1 = scaffold site
        --> 0 = not a scaffold site

        Arguments:
            coords - list of coordinates representing every lattice
            site
        """
        lattice = deepcopy(coords)
        y_max = lattice[:, 1].max()
        x_max = lattice[:, 0].max()

        grid = np.zeros((y_max + 1, x_max + 1)).astype(int)
        for point in lattice:
            grid[point[1], point[0]] = 1
        return grid

    @staticmethod
    def quantise_rows(array, poss_cross, padding):  # Returns array
        """
        Ensure first and last scaffold site in the row correlate to
        a position where the half turn crossover most likely occurs.
        Values are rounded up OR down to the closest half turn location.

        Arguments:
            array - binary array, where 1's represent lattice sites

        i.e. a row of length 9 will round down to 5 lattice sites
        or 43 will round up to 46
        """
        grid = deepcopy(array)

        # Add a border of 16 0's around the lattice
        grid = np.pad(grid, pad_width=padding, mode="constant", constant_values=0)

        # Find the number of nucleotide sites per row and store in a numpy array
        nt_per_row = []
        for row in range(len(grid)):
            nt_per_row.append(np.sum(grid[row]))
        nt_per_row = np.array(nt_per_row)

        # Find closest crossover and modify rows ensure first/last sites are crossover sites

        # values 1-4 round up
        nt_per_row_round_1 = [5 if 0 < i < 5 else i for i in nt_per_row]
        # other values round to closest half turn
        closest_crossover = lambda x: find_closest(poss_cross, x)
        nt_per_row_round_2 = np.array(
            list(map(closest_crossover, nt_per_row_round_1)))
        nt_per_row_diff = nt_per_row_round_2 - nt_per_row

        grid = modify_lattice_row(grid, nt_per_row_diff)

        return grid

    @staticmethod
    def straighten_edges(array, sf: int):
        """
        Algorithm to straighten out edges of the lattice
        
        Arguments:
            array - binary array of 1's and 0's
            sf - straightening factor

        Method:
        Shift each row to begin at a multiple of sf (straightening_factor)

        """
        grid = deepcopy(array)
        for row in range(np.shape(grid)[0]):
            if np.sum(grid[row]) == 0:  # to account for padding
                continue
            else:
                starting_location = int(np.argwhere(grid[row])[0])
                rounded_location = round_to_multiple(starting_location, sf, 0)
                change_in_location = rounded_location - starting_location
                grid[row] = np.roll(grid[row], change_in_location)
        return grid

    @staticmethod
    def connect_rows(array, start_side, half_turn_locations, bp_per_turn):
        """
        Algorithm to ensure every row connects to the next

        Arguments:
            array - binary array of 1's and 0's
        """

        grid = deepcopy(array)
        # half turn nt sizes [5,16,26,37,47, etc] instead of half turn indexes [0,4,15, etc]
        poss_cross = np.add(half_turn_locations[1:], 1)
        width_tracker = []

        for index in range(np.shape(grid)[0] - 3):

            row0 = grid[index]
            row1 = grid[index+1]
            row2 = grid[index+2]
            row3 = grid[index+3]
            row_width = [np.sum(row0),np.sum(row1), np.sum(row2), np.sum(row3)]
            width_tracker.append(int(row_width[0]))
            if 0 in row_width[0:2]: continue 

            # R represents the side where the first crossover occurs
            # Hence, if self.start_side = "left"
            # for 0,2,4 etc -> R = True, otherwise R = None
            if start_side == "left":
                R = True if index % 2 == 0 else None
            else:
                R = None if index % 2 == 0 else True

            # calculate current difference between end/start points of row0 & row1
            if R:
                right0 = int(np.argwhere(row0)[-1])
                right1 = int(np.argwhere(row1)[-1])                
                if right0 - right1 == 0: continue
                
                extra_turns_right = int(round(abs(right1-right0) / bp_per_turn, 0))
            else:
                left0 = int(np.argwhere(row0)[0])
                left1 = int(np.argwhere(row1)[0])
                if left0 - left1 == 0: continue

                extra_turns_left = int(round(abs(left1-left0) / bp_per_turn, 0))

            ### Initialise some useful values
            # find index in poss_cross correlating to no. of possible crossovers in that row
            cross_idx_bottom = bisect_left(poss_cross, row_width[0])
            cross_idx_top = bisect_left(poss_cross, row_width[1])
            TminusB = row_width[1] - row_width[0]
            top_bigger = True if TminusB > 0 else False
            extra_turns = cross_idx_bottom - cross_idx_top  # in whole row
            row0_diff = 0
            row1_diff = 0

            # Same size rows, just shift them left or right
            if TminusB == 0:
                if R:
                    row1_roll = right0 - right1
                else:
                    row1_roll = left0 - left1
                row1 = np.roll(row1, row1_roll)

            elif abs(extra_turns) > 1:
                if not top_bigger:
                    # shorten the end of row 0
                    if row_width[0] <= poss_cross[3] and not all_same(width_tracker[-3:-1]):
                        row0_diff = poss_cross[cross_idx_bottom - 1] - row_width[0]

                    else:
                        if R:
                            row0_diff = poss_cross[cross_idx_bottom - extra_turns_right] - row_width[0] 
                        else:
                            row0_diff = poss_cross[cross_idx_bottom - extra_turns_left] - row_width[0]
                else:
                    # lengthen the end of row 0
                    if row_width[0] <= poss_cross[3] and not all_same(width_tracker[-3:-1]):
                        row0_diff = poss_cross[cross_idx_bottom + 1] - row_width[0]
                    else:
                        if R:
                            row0_diff = poss_cross[cross_idx_bottom + extra_turns_right] - row_width[0] 
                        else:
                            row0_diff = poss_cross[cross_idx_bottom + extra_turns_left] - row_width[0]

                row0 = modify_lattice_row(row0, row0_diff, side(0, R))
                row1 = modify_lattice_row(row1, row1_diff, side(1, R))

                if R:
                    row1_roll = right0 - right1 + row0_diff
                    row1 = np.roll(row1, row1_roll)
                elif not R:
                    row1_roll = left0 - left1 - row0_diff
                    row1 = np.roll(row1, row1_roll)

            elif abs(extra_turns) == 1:
                if cross_idx_bottom >= 1:  #16 or bigger
                    if top_bigger:
                        #if row1 is smaller than the last few rows and the next few rows are equal
                        # if not row_width[1] >= max(width_tracker[-5:-1]) and row_width[1] == row_width[2] == row_width[3]:
                        #     # remove 1 to row1 -> so bigger shapes are less skewed
                        #     row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                        #     row1 = modify_lattice_row(row1, row1_diff, side(1,R)) 
                        
                        # if the next 2 rows are equal and 
                        if row_width[1] == row_width[2] and not row_width[1] > max(width_tracker[-10:-1]): 
                            row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                            row1 = modify_lattice_row(row1, row1_diff, side(1,R)) 
                            # pass
                        
                        elif row_width[2] in [row_width[1], 0]:
                            if max(width_tracker[-10:-1]) >= row_width[0]:
                                # shorten row1 by 1 crossover
                                row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                                row1 = modify_lattice_row(row1, row1_diff, side(1, R))

                    elif not top_bigger:
                        # if the next 2 rows and the previous row was the same and the current row isn't bigger than the last few rows
                        if all_same([width_tracker[-2], row_width[1], row_width[2]]):
                            #if row_width[0] > max(width_tracker[0:-1]):
                            # shorten end of row 0 by 1 crossover
                            row0_diff = poss_cross[cross_idx_bottom - 1] - row_width[0]
                            row0 = modify_lattice_row(row0, row0_diff, side(0, R))

                        # where the last 2 rows have been equal and the next 2 rows are equal (or equal to 0)
                        elif min(width_tracker[-3:-1]) in [max(width_tracker[-3:-1]), 0] and (all_same(row_width[1:3]) or row_width[3] == row_width[2] == 0):
                            if cross_idx_top > 0:
                                # shorten ends of row 0/1 by 1 crossover
                                row0_diff = poss_cross[cross_idx_bottom - 1] - row_width[0]
                                row0 = modify_lattice_row(row0, row0_diff, side(0, R))
                                row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                                row1 = modify_lattice_row(row1, row1_diff, side(1, R))
                            elif cross_idx_top < 1:
                                # shorten ends of row 0/1 by 1 crossover
                                row0_diff = poss_cross[cross_idx_bottom - 1] - row_width[0]
                                row0 = modify_lattice_row(row0, row0_diff, side(0, R))

                        # when row0 is size ~36 and the next 2 rows are equal, lengthen row1 by 1
                        elif cross_idx_bottom == 3 and row_width[1] == row_width[2]:
                            # lengthen end of row1 by 1 crossover
                            row1_diff = poss_cross[cross_idx_top + 1] - row_width[1]
                            row1 = modify_lattice_row(row1, row1_diff, side(1, R))

                        # where the row length is > 45 and the next row isnt bigger than the last 2 rows
                        # add 1 to row1 -> so bigger shapes are less skewed
                        # elif cross_idx_bottom > 4 and not row_width[1] > max(width_tracker[-3:-1]):
                        #     row1_diff = poss_cross[cross_idx_top + 1] - row_width[1]
                        #     row1 = modify_lattice_row(row1, row1_diff, side(1, R))

                elif cross_idx_bottom == 0:  #5
                    if top_bigger:
                        if row_width[2] < row_width[1] or row_width[2] <= poss_cross[1]:
                            # row 1: 16 --> 5
                            row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                            row1 = modify_lattice_row(row1, row1_diff, side(1, R))
                    elif not top_bigger:
                        #just shift
                        pass

                if R:
                    row1_roll = right0 - right1 + row0_diff
                    row1 = np.roll(row1, row1_roll)
                elif not R:
                    row1_roll = left0 - left1 - row0_diff
                    row1 = np.roll(row1, row1_roll)

            # UPDATE ROWS IN GRID
            grid[index] = row0
            grid[index + 1] = row1
            width_tracker[index] = np.sum(row0)

        return grid

    @staticmethod
    def array_to_coords(array):
        """Returns lattice points as a set of coordinates"""
        # Remove padding and find coordinates of updated system
        grid = deepcopy(array)
        y_min = np.argwhere(grid)[:, 0].min()
        y_max = np.argwhere(grid)[:, 0].max()
        x_min = np.argwhere(grid)[:, 1].min()
        x_max = np.argwhere(grid)[:, 1].max()

        grid = grid[y_min:y_max + 1, x_min:x_max + 1]
        coords = np.argwhere(grid)
        # swap columns, np.argwhere returns x and y the wrong way around
        coords[:, 0], coords[:, 1] = coords[:, 1], coords[:, 0].copy()
        return coords

    @staticmethod
    def get_crossovers(array, start_side, poss_cross):
        """ Returns a binary array, where 1 represents a potential crossover location"""
        grid = deepcopy(array)
        # copy the size of grid but fill it with zeros
        crossovers_array = np.zeros(np.shape(grid))

        sides = {"left": "right", "right": "left"}
        side = start_side  # begin at the left / rightside
        for row in range(np.shape(grid)[0]):
            # to account for padding, calc no. of lattice sites on row
            lattice_sites = np.sum(grid[row])
            if lattice_sites == 0:  # i.e. when row is just padding
                continue

            # find x of first/last lattice site in the row
            left_bound = np.argwhere(grid[row])[0]
            right_bound = np.argwhere(grid[row])[-1]

            for bp in poss_cross:
                if bp > lattice_sites:
                    continue
                if side == "left":
                    crossovers_array[row, left_bound + bp] = (1 * grid[row, left_bound + bp])
                else:  # if side == right
                    crossovers_array[row, right_bound - bp] = (1 * grid[row, right_bound - bp])

            side = sides[side]

        return crossovers_array

    def plotPolygon(self, ax, nodes: np.ndarray, coords: bool):
        polygon = deepcopy(self.polygon_array)
        if np.shape(nodes)[1] not in [2,3]: # if array
            nodes = self.array_to_coords(nodes)
        
        x_min = polygon[:, 0].min()
        y_min = polygon[:, 1].min()

        # shift polygon to (0,0)
        polygon[:, 0] -= x_min
        polygon[:, 1] -= y_min

        # Normalise values to integer values, e.g. 0.34 -> 1
        polygon[:, 0] /= self.x_spacing
        polygon[:, 1] /= self.y_spacing

        # add first value to start
        polygon = np.vstack((polygon, polygon[0]))

        if coords:
            # Center polygon
            polygon[:, 0] += (nodes[:, 0].max() - polygon[:, 0].max()) / 2
        else:
            # flip
            polygon[:, 1] = polygon[:, 1].max() - polygon[:, 1]
            # add padding to match
            polygon += [self.padding, self.padding, 0]

        ax.plot(polygon[:, 0], polygon[:, 1], "r--", linewidth = 1, alpha = 0.5)

        return ax

    def plot(self,
             plot_these: List[np.ndarray],
             ax: plt.Axes = None,
             fout: str = None,
             poly = True,
             root: str = "",
             title: str = None,
             aspect: int = 3,
             ticks: list = [10,5],
             grid: bool = False):
        """
        Plot the lattice coordinates or lattice as an array

        Arguments:
            plot_these - will take up to 2 sets of lattice and crossover coordinates
            ax - an axes
            fout - saves the plot as a png with the name of `fout`
            poly - True will plot the initial polygon vertices
            root - directory to save file in, defaults to saving in current directory
            title - this is the title of the plot
            aspect - this is the aspect ratio of the x and y axes

        """
        assert len(plot_these) <= 4, "Max no. of plots on one axis reached, 4 or less"

        if not ax:
            fig, ax = plt.subplots()
        if grid:
            plt.grid(True)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(4)
        ax.xaxis.set_major_locator(MultipleLocator(ticks[0]))
        ax.yaxis.set_major_locator(MultipleLocator(ticks[1]))
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("No. of strands")

        point_style = itertools.cycle(["ko", "b.","r.","cP"])
        point_size = itertools.cycle([0.5, 2.5])
        for points in plot_these:
            if np.shape(points)[1] not in [2,3]: # if array
                nodes = self.array_to_coords(points)
            else:
                nodes = points
            # Lattice sites then crossover sites
            ax.plot(nodes[:, 0], nodes[:, 1], next(point_style), ms=next(point_size), alpha=0.25)
        
        if poly:
            self.plotPolygon(ax, plot_these[0], coords=True)
        if title:
            ax.set_title(f"{title}")

        plt.gca().set_aspect(aspect)
        if fout:
            plt.savefig(f"{root}{fout}.png", dpi=500)
        if not ax:
            plt.show()

class DNASnake(Lattice):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        """
        Subclass of Lattice. Inherits all atributes and methods from this parent class.

        When a DNASnake object is generated, it refers to a DNA Scaffold 
        with a maximum of 2 crossovers per "row" in the lattice
        """
        self.quantised_array = self.quantise_rows(self.intersected_array, self.poss_cross, self.padding)
        # self.straightened_array = self.straighten_edges(self.quantised_array, self.straightening_factor)
        self.connected_array = self.connect_rows(self.quantised_array, 
                                                self.start_side, 
                                                self.poss_cross, 
                                                self.bp_per_turn)

        # Final Lattice and crossovers
        self.final_array = self.connected_array
        self.crossover_array = self.get_crossovers(self.final_array, self.start_side, self.poss_cross)

        # Lattice and Crossovers as coordinates
        self.final_coords = self.array_to_coords(self.final_array)
        self.crossover_coords = self.array_to_coords(self.crossover_array)

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
