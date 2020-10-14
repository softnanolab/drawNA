import numpy as np

from drawNA.tools import DNANode, DNAEdge

def test_Node():
    node = DNANode(np.array([0., 1., 0.]))
    assert node.angle == None
    node.vector_3p = np.array([0., 1., 0.])
    node.vector_5p = np.array([1., 0., 0.])
    assert node.angle == np.pi / 2
    return

def test_Edge():
    nodes = [
        DNANode(np.array([0., 1., 0.])),
        DNANode(np.array([0., 3., 0.]))
    ]
    edge = DNAEdge(*nodes)
    assert edge.length == 2.0
    return

def test_complex():
    nodes = [
        DNANode(np.array([0., 1., 0.])),
        DNANode(np.array([0., 3., 0.])),
        DNANode(np.array([3., 3., 0.])),
        DNANode(np.array([3., 1., 0.])),
    ]

    edges = [DNAEdge(nodes[i-1], node) for i, node in enumerate(nodes)]
    for i, node in enumerate(nodes):
        print(f"{i}: {node.vector_3p}, {node.vector_5p}")
        assert node.angle == np.pi / 2
    return

if __name__ == '__main__':
    test_Node()
    test_Edge()
    test_complex()