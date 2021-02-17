from drawNA.lattice.route import LatticeRoute
from drawNA.lattice import LatticeRoute
from drawNA.polygons import BoundaryPolygon
from drawNA.oxdna import System
from ...sandbox.staples_from_route import StapleBaseClass, StapleContainer, StaplingAlgorithm1

import numpy as np

def test_run_staple_algorithm_1(domain_size = 15) -> StapleBaseClass:
    # Create Scaffold strand object (i.e. LatticeRoute)
    polygon = BoundaryPolygon(polygon_vertices) # this step can be omitted in theory but shouldn't be
    scaffold = polygon.dna_snake(start_side="left", grid_size = [0.34, 2])
    route = scaffold.route()

    # Generate staples for all nucleotides on scaffold
    print("Test: Generating staples for all nucleotides on the scaffold")
    stapled_scaffold = StaplingAlgorithm1(route, domain_size = domain_size)
    return stapled_scaffold

def test_create_staple_container(stapled_scaffold: StapleBaseClass) -> StapleContainer:
    # Generate oxDNA strand objects for every staple and add to a container
    print("Test: Adding staples to container")
    strand_container = stapled_scaffold.generate_container()
    return strand_container

def test_return_system(container: StapleContainer) -> System:
    print("Test: Adding staples to an oxDNA system...")
    return container.system()
    


if __name__ == "__main__":
    # Create Scaffold strand object (i.e. LatticeRoute)
    polygon_vertices = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]]) # trapezium
    polygon = BoundaryPolygon(polygon_vertices) # this step can be omitted in theory but shouldn't be
    scaffold = polygon.dna_snake(start_side="left", grid_size = [0.34, 2])
    route = scaffold.route()

    # Create Stapled Scaffold Object
    stapled_scaffold = test_run_staple_algorithm_1(route)
    # Add all strands into a StapleContainer
    strand_container = test_create_staple_container(stapled_scaffold)
    # Return oxDNA system
    system = test_return_system(strand_container)
    system.write_oxDNA(prefix = "test_trapezium")
