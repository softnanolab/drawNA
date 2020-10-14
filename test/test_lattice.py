from os import path

import numpy as np

from drawNA.lattice import LatticeNode, LatticeEdge, LatticeRoute
from drawNA.oxdna import System

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])

def test_node():
    node = LatticeNode(np.array([0., 1., 0.]))
    assert node.angle == None
    node.vector_3p = np.array([0., 1., 0.])
    node.vector_5p = np.array([1., 0., 0.])
    assert node.angle == np.pi / 2
    return

def test_edge():
    nodes = [
        LatticeNode(np.array([0., 1., 0.])),
        LatticeNode(np.array([0., 3., 0.]))
    ]
    edge = LatticeEdge(*nodes)
    assert edge.length == 2.0
    return

def test_route():
    nodes = [
        LatticeNode(np.array([0., 10., 0.])),
        LatticeNode(np.array([0., 30., 0.])),
        LatticeNode(np.array([30., 30., 0.])),
        LatticeNode(np.array([30., 10., 0.])),
    ]
    route = LatticeRoute(nodes)
    assert len(route.edges) == 3
    system = System(np.array([50., 50., 50.]))
    system.add_strand(route)
    assert len(system.strands) == 1
    system = route.system(box=np.array([50., 50., 50.]))
    assert len(system.strands) == 1
    system.write_oxDNA(root=ROOT)
    return route

if __name__=='__main__':
    test_node()
    test_edge()
    route = test_route()
    system = route.system()
    
    #print(route)
    #route.plot()