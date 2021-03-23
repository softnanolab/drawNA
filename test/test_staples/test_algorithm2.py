from drawNA.lattice.route import LatticeRoute
from drawNA.lattice import LatticeRoute
from drawNA.polygons import BoundaryPolygon
from drawNA.oxdna import System
from sandbox.staples_from_route import ConfGenSplitDoubleDomains, StapleBaseClass, StaplingAlgorithm2, ConfGenSplitDoubleDomains

import numpy as np
import pytest
from os import path

### ---
# This is an algorithm specifically for shapes with aligned left and right sides, i.e. square
### ---

# 1. Arrange
@pytest.fixture
def scaffold() -> LatticeRoute: 
    # Create Scaffold strand object (i.e. LatticeRoute)
    square = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])*np.array([3.5,2.5,1])
    rectangle = square*[8,5,1]
    polygon = BoundaryPolygon(rectangle) # this step can be omitted in theory but shouldn't be
    scaffold = polygon.dna_snake(start_side="left", grid_size = [0.34, 2])
    route = scaffold.route()
    return route

# 2. Act
@pytest.fixture
def stapled_scaffold(scaffold: LatticeRoute) -> StapleBaseClass:    
    # Generate staples for all nucleotides on scaffold
    stapled_scaffold = StaplingAlgorithm2(scaffold)
    return stapled_scaffold

# 3. Assert
def test_container(stapled_scaffold: StapleBaseClass):
    # Generate oxDNA strand objects for every staple and add to a container
    container = stapled_scaffold.generate_origami()
    # Generate oxDNA system filled with strand objects
    system = container.system()
    system.write_oxDNA(prefix = "test_square_algo_2")

    # Assertion statements
    assert container.n_staples == stapled_scaffold.staple_ID - 1
    ## add lots more later

    return container

def test_configuration(stapled_scaffold: StapleBaseClass):
    # Generate oxDNA strand objects for every staple and add to a container
    container = stapled_scaffold.generate_origami()
    # Generate two single domain strands for all double domain strands which are not at the edges of the scaffold
    generator = ConfGenSplitDoubleDomains(staple_strands = container.staples, staple_base_class = stapled_scaffold)
    
    # First configuration is the original configuration, hence n_staples of all others should be one greater than this
    n_staples_in_each_config = np.array([conf.n_staples for conf in generator.configurations])
    assert np.all(n_staples_in_each_config[1:] == n_staples_in_each_config[0] + 1 )
    
    # 8 inner staples, 6 of which are double-domained. Therefore, expect original + 6 configurations.
    assert len(generator.configurations) == 7
    
    # Test Write To File Function
    generator.write_to_file(name = "batch1") 


## Test in jupyter interactive terminal
# !pytest test_algorithm_1.py

# polygon_vertices = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]])*6 # trapezium
# polygon = BoundaryPolygon(polygon_vertices) # this step can be omitted in theory but shouldn't be
# scaffold = polygon.dna_snake(start_side="left", grid_size = [0.34, 2])
# route = scaffold.route()
# route.plot()
# stapled_scaf = StaplingAlgorithm1(route, domain_size = 15)
# container = stapled_scaf.generate_container()