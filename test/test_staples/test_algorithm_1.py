from drawNA.lattice.route import LatticeRoute
from drawNA.lattice import LatticeRoute
from drawNA.polygons import BoundaryPolygon
from drawNA.oxdna import System
from sandbox.staples_from_route import StapleBaseClass, StapleContainer, StaplingAlgorithm1

import numpy as np
import pytest

# 1. Arrange
@pytest.fixture
def scaffold() -> LatticeRoute: 
    # Create Scaffold strand object (i.e. LatticeRoute)
    polygon_vertices = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]])*6 # trapezium
    polygon = BoundaryPolygon(polygon_vertices) # this step can be omitted in theory but shouldn't be
    scaffold = polygon.dna_snake(start_side="left", grid_size = [0.34, 2])
    route = scaffold.route()
    route.plot()
    return route

# 2. Act
@pytest.fixture
def stapled_scaffold(scaffold: LatticeRoute) -> StapleBaseClass:    
    # Generate staples for all nucleotides on scaffold
    stapled_scaffold = StaplingAlgorithm1(scaffold, domain_size = 15)
    return stapled_scaffold

# 3. Assert
def test_container(stapled_scaffold: StapleBaseClass):
    # Generate oxDNA strand objects for every staple and add to a container
    container = stapled_scaffold.generate_container()

    # Generate oxDNA system filled with strand objects
    system = container.system()
    system.write_oxDNA(prefix = "test_trapezium")

    # Assertion statements
    assert container.n_staples == stapled_scaffold.staple_ID - 1
    ## add lots more later

    return container


## Test in jupyter interactive terminal
# !pytest test_algorithm_1.py

# polygon_vertices = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]])*6 # trapezium
# polygon = BoundaryPolygon(polygon_vertices) # this step can be omitted in theory but shouldn't be
# scaffold = polygon.dna_snake(start_side="left", grid_size = [0.34, 2])
# route = scaffold.route()
# route.plot()
# stapled_scaf = StaplingAlgorithm1(route, domain_size = 15)
# container = stapled_scaf.generate_container()