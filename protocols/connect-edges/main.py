#!/usr/bin/env python
"""Manual creation of DNAEdge instances that are joined correctly"""
import numpy as np
from drawNA.lattice import LatticeEdge, LatticeNode
from drawNA.oxdna import System, Strand
import drawNA.oxdna.utils as ox
from drawNA.oxdna.strand import generate_helix

def with_edges():
    nodes = [
        LatticeNode(np.array([0., 0., 0.])),
        LatticeNode(np.array([0., 5., 0.])),
        LatticeNode(np.array([0., 10., 0.])),
    ]
    edges = [
        LatticeEdge(nodes[0], nodes[1]),
        LatticeEdge(nodes[1], nodes[2], theta=np.pi),
    ]
    
    for i, edge in enumerate(edges):
        print(f"Edge[{i}]: {edge.summary}")
        print(f"\t0x7f...{hex(id(edge.vertices[0]))[-6:]}, 0x7f...{hex(id(edge.vertices[1]))[-6:]}")
        strands = edge.strand(double=True)
        print(f"\t{strands[0].nucleotides[0].pos_com}")
        print(f"\t{strands[0].nucleotides[-1].make_5p(None).pos_com}")
    return strands

def with_strands(N : int = 10):
    a1 = np.array([1., 0., 0.])
    start = np.array([0., 0., 0.]) + a1 * 0.6
    direction = np.array([0., 1., 0.])
    strands = generate_helix(
        n=N,
        start_position=start,
        direction=direction,
        a1=a1,
        initial_rotation=0.0,
    )
    new_position = start + (N) * 0.39 * direction
    new_direction = np.array([0.0, 1.0, 1.0]) / 2 ** 0.5
    new_a1 = np.array([0.0, -1.0, 1.0]) / 2 ** 0.5
    #print(np.dot(new_a1, new_direction))
    rotation = -0.06
    #rotation = 0.0

    strands += generate_helix(
        n=N,
        start_position=new_position,
        direction=new_direction,
        a1=new_a1,
        initial_rotation=rotation,
    )
    strand = Strand(strands[0].nucleotides + strands[1].nucleotides)
    return [strand]

def main(strands : bool = True, edges : bool = False):
    system = System(np.array([30., 30., 30.]))
    if strands:
        system.add_strands(with_strands())
    if edges:
        system.add_strands(with_edges())
    system.write_oxDNA('connect')
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    main(**vars(parser.parse_args()))