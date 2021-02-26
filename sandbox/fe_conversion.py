from typing import List

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

def main():
    system = Converter([50, 50, 50])
    strands = generate_helix(21, double=True)
    system.add_strands(strands)
    print(system.dataframe)
    nodes = system.create_nodes()
    print(nodes)
    strands = generate_helix(5, double=True)
    system = System([50, 50, 50])
    system.add_strands(strands)
    system.nucleotides
    print(json.dumps({hex(id(i)): f"{i.__repr__()}: [{i.across}|{i._across.index}] {hex(id(i._across))}" for i in system.nucleotides}, indent=2))

    
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

if __name__ == '__main__':
    import json
    main()
    