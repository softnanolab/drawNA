#!/usr/bin/env python
"""
Protocols for demonstrating the usage of the LatticeRoute class.
"""
import numpy as np

from drawNA.lattice import LatticeEdge, LatticeNode, LatticeRoute
from drawNA.oxdna import Nucleotide, Strand, System
from drawNA.oxdna.strand import generate_helix

def long_strand() -> Strand:
    strands = generate_helix(20)
    return strands[0]

def long_route() -> LatticeRoute:
    nodes = [
        LatticeNode(np.array([0., 0., 10.])),
        LatticeNode(np.array([0, 20. * 0.42, 10.]))
    ]
    route = LatticeRoute(nodes)
    return route

def multi_route() -> LatticeRoute:
    nodes = [
        LatticeNode(np.array([0., 0., 20.])),
        LatticeNode(np.array([0., 10. * 0.42, 20.])),
        LatticeNode(np.array([0., 20. * 0.42, 20.]))
    ]
    route = LatticeRoute(nodes)
    return route

def strand_length(strand : Strand) -> float:
    return np.linalg.norm(strand.nucleotides[-1].pos_base - strand.nucleotides[0].pos_base)

def main():
    strand = long_strand()
    print(f"Strand:\n    {strand}")
    print(f"    Length between bases of the end nucleotides: {strand_length(strand)}")
    route = long_route()
    print(f"\nSingle-Edge Route:\n    {route}")
    print(f"    Length between bases of the end nucleotides: {strand_length(route)}")
    multi = multi_route()
    print(f"\nDouble-Edge Route:\n    {multi}")
    print(f"    Length between bases of the end nucleotides: {strand_length(multi)}")
    system = System(np.array([50., 50., 50.]))
    system.add_strands([strand, route, multi])
    system.write_oxDNA('route')
    return

if __name__ == '__main__':
    print('Running easy-route...\n')
    main()