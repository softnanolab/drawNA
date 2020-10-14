import numpy as np

from ..tools import DNANode

class LatticeNode(DNANode):
    def __init__(self, position : np.ndarray):
        super().__init__(position)