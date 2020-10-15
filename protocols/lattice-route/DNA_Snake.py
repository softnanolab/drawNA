from drawNA.oxdna import System
from drawNA.lattice import LatticeRoute, DNASnake
from drawNA.polygons import BoundaryPolygon

import numpy as np

from os import path

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])

def generate(polygon_vertices: np.ndarray, title: str = "polygon_name", DNAout: bool = False, PLOTout: bool = False):
    print(f"{title}: Making polygon...")
    polygon = BoundaryPolygon(polygon_vertices)
    print(f"{title}: ...constructing scaffolding lattice...")
    lattice = polygon.dna_snake(straightening_factor=5, start_side="left", grid_size = [0.34, 2.5])
    print(f"{title}: ...calculating route.")
    route = lattice.route()

    
    if PLOTout:
        print(f"{title}: Generating and saving lattice plot.")
        plot_list = [lattice.quantised_array, lattice.crossover_coords]
        lattice.plot(plot_list, poly = True, fout = title, root = ROOT + "/lattice-plot_out/", aspect = 6)
    route.plot()

    if DNAout:
        system = route.system()
        print(f"{title}: Generating and saving oxDNA files.")
        system.write_oxDNA(prefix = title, root=ROOT+"/oxDNA_out/")
        print(f"{title}: Generating and saving LAMMPS files.")
        system.write_LAMMPS(prefix = title, root=ROOT+"/LAMMPS_out/")


def main():
    """ Polygon Vertices """
    square = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]])
    trap = np.array([[0.,0.,0.],[1.5,6.,0.],[8.5,6.,0.],[10.,0.,0.]])
    trapREV = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]])
    no = 0.85
    hexagon = np.array(
        [[0, 1, 0], [no, 0.5, 0], [no, -0.5, 0], [0, -1, 0], [-no, -0.5, 0], [-no, 0.5, 0]])  
    plus = np.array([
        [1.,0.,0.], [2.,0.,0.], [2.,1.,0.], [3.,1.,0.], [3.,2.,0.], [2.,2.,0.],
        [2.,3.,0.], [1.,3.,0.], [1.,2.,0.], [0.,2.,0.], [0.,1.,0.], [1.,1.,0.]])
    diamond = np.array([[1.,0.,0.],[2.,1.,0.],[1.,2.,0.],[0.,1.,0.]])
    triangle = np.array([[0,0,0],[5,9,0],[10,0,0]])
    stacked_I = np.array([
    [0.,0.,0.],[3.,0.,0.],[3.,1.,0.],[2.,1.,0.], [2.,2.,0.],[3.,2.,0.],
    [3.,3.,0.],[2.,3.,0.],[2.,4.,0.],[3.,4.,0.],[3.,5.,0.],[0.,5.,0.],[0.,4.,0.],[1.,4.,0.],
    [1.,3.,0.],[0.,3.,0.], [0.,2.,0.],[1.,2.,0.],[1.,1.,0.],[0.,1.,0.]
    ])
    octagon = np.array([[1.0,-0.85,0],[2.0, -0.85, 0],[2.5, 0, 0],[2.5, 1,0],
    [2.0, 1.85,0],[1.0,1.85,0],[0.5,1,0],[0.5, 0.0,0.0]])

    """ Generating Scaffold Strands """
    print("\n ----- Generating scaffold for a square: \n")
    generate(square*4, title = "square", DNAout=True, PLOTout=True)

    print("\n ----- Generating scaffold for a trapezium: \n")
    generate(trap*6, title="trapezium", DNAout=True, PLOTout=True)

    print("\n ----- Generating scaffold for a hexagon: \n")
    generate(hexagon*15, title="hexagon", DNAout=True, PLOTout=True)

    print("\n ----- Generating scaffold for a plus: \n")
    generate(plus*10, title="plus", DNAout=True, PLOTout=True)

    print("\n ----- Generating scaffold for a modular shape 'stacked_I': \n")    
    generate(stacked_I*12, title="stacked_I", DNAout=True, PLOTout=True)

    """ Uncomment/add/modify generate commands to try other shapes """
    # generate(diamond*10, title="diamond", DNAout=True, PLOTout=True)
    # generate(trapREV*5, title="trapezium_rev",DNAout=True, PLOTout=True)
    # generate(octagon*12, title="oct", DNAout=True, PLOTout=True)
    # generate(triangle*2, title="triangle",DNAout=True, PLOTout=True)

if __name__ == "__main__":
    main()
    print("\n NOTE: See DNA_Snake.py to understand how to generate a scaffold for a given shape")