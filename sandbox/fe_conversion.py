from typing import List

from drawNA.oxdna import Nucleotide, Strand, System
from drawNA.oxdna import generate_helix

from Spring_system import Node, Spring, Beam, Connector # type: ignore

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

def main():
    system = Converter([50, 50, 50])
    strands = generate_helix(5, double=True)
    system.add_strands(strands)
    print(system.dataframe)
    nodes = system.create_nodes()
    print(nodes)
    return

if __name__ == '__main__':
    import json
    main()
    strands = generate_helix(5, double=True)
    system = System([50, 50, 50])
    system.add_strands(strands)
    system.nucleotides
    print(json.dumps({hex(id(i)): f"{i.__repr__()}: [{i.across}|{i._across.index}] {hex(id(i._across))}" for i in system.nucleotides}, indent=2))