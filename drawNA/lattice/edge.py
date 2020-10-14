from ..tools import DNAEdge

class LatticeEdge(DNAEdge):
    def __init__(self, vertex_1, vertex_2, **kwargs):
        super().__init__(vertex_1, vertex_2, **kwargs)