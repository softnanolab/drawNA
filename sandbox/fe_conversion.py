import json

from typing import List

import pandas as pd
import matplotlib.pyplot as plt

from drawNA.oxdna import Nucleotide, Strand, System
from drawNA.oxdna import generate_helix

from Spring_system import Node, Spring, Beam, Connector # type: ignore


from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class NucleotideSingle(Node):
    def __init__(self, nt_1: Nucleotide):
        self.nt_1 = nt_1
        self.nt_2 = None
        super().__init__()
        self.set_position(
            nt_1.pos_com, # e0
            nt_1._a1, # e1
            nt_1._a2, # e2
            nt_1._a3, # e3
        )

class NucleotidePair(Node):
    def __init__(self, nt_1: Nucleotide, nt_2: Nucleotide):
        self.nt_1 = nt_1
        self.nt_2 = nt_2
        super().__init__()
        self.set_position(
            nt_1.pos_com + 0.5 * (nt_2.pos_com - nt_1.pos_com), # e0
            nt_1._a1, # e1
            nt_1._a2, # e2
            nt_1._a3, # e3
        )

class Converter(System):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def create_nodes(self) -> List[Node]:
        result = [] # list of Node instances
        completed = [] # list of completed indices

        nucleotides = self.nucleotides

        for nt in nucleotides:
            if nt.index in completed: continue
            if nt._across:
                node = NucleotidePair(nt, nt._across)
                result.append(node)
                completed += [nt.index, nt.across]

            else:
                node = NucleotideSingle(nt)
                result.append(node)
                completed.append(nt.index)

        return result

    def create_connectors(self) -> List[Connector]:
        connectors = []
        # use self.bonds property to iterate over all bonds
        # each bond should be assigned a connector

        # 
        return

    @property
    def system_FE(self) -> list:
        return self._nodes + self._connectors

def add_node(ax: plt.Axes, node: Node):
    position = node.position[0]
    e1 = node.position[1]
    e2 = node.position[2]
    e3 = node.position[3]
    ax.plot([position[0]], [position[1]], [position[2]], 'ko', markersize=10)
    
    colours = 'r', 'g', 'b'
    for i, v in enumerate([e1, e2, e3]):
        print(i, v)
        vector = 0.05 * v
        ax.quiver(*position, *vector, color=colours[i], length=2)
        continue
        arrow = Arrow3D(
            [position[0], position[0] + vector[0]],
            [position[1], position[1] + vector[1]],
            [position[2], position[2] + vector[2]],
            arrowstyle="-|>", color=colours[i], lw=5, mutation_scale=20,
        )
        ax.add_artist(arrow)
    return

def create_converter(kind: str) -> Converter:
    def _duplex() -> List[Strand]:
        return generate_helix(5, 'ATCGG', double=True)

    def _nick() -> List[Strand]:
        strands = generate_helix(6, 'AATTCC', double=True)
        main_strand = strands[0]
        single_strands = [
            Strand(
                nucleotides=strands[1].nucleotides[i*3: (i+1)*3]
            ) for i in range(2)
        ]
        return [main_strand] + single_strands

    def _ssDNA() -> List[Strand]:
        strands = generate_helix(11, 'GTAGTAATGGG', double=True)
        main_strand = strands[0]
        for nt in main_strand.nucleotides[3:8]:
            nt.across = None
        single_strands = [
            Strand(
                nucleotides=strands[1].nucleotides[i: i+3]
            ) for i in (0, 8)
        ]
        return [main_strand] + single_strands

    def _double_crossover() -> List[Strand]:
        return

    def _single_crossover() -> List[Strand]:
        return

    def _open_nick() -> List[Strand]:
        return

    def _bulge() -> List[Strand]:
        return

    _dispatch = {
        'duplex': _duplex,
        'nick': _nick,
        'ssDNA': _ssDNA,
        'double_crossover': _double_crossover,
        'single_crossover': _single_crossover,
        'open_nick': _open_nick,
        'bulge': _bulge,
    }

    system = Converter([20, 20, 20])
    strands = _dispatch[kind]()
    if strands:
        system.add_strands(strands)
    return system

def visualise_nodes(nodes: List[Node]):
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.autoscale(enable=True, axis='both', tight=True)
    
    ax.set_xlim((0, 2))
    ax.set_ylim((0, 2))
    ax.set_zlim((0, 2))
    for node in nodes:
        add_node(ax, node)

    plt.tight_layout()
    plt.show()
    return

def main():
    names = [
        'duplex',
        'nick',
        'ssDNA',
        'double_crossover',
        'single_crossover',
        'open_nick',
        'bulge',
    ]

    # create a dictionary of name with its corresponding system
    systems = {name: create_converter(name) for name in names}
    return systems

if __name__ == '__main__':
    systems = main()
    print(systems['ssDNA'].dataframe)
    