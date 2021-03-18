import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import itertools
from typing import Tuple
import pandas as pd
'''
Node class creates a node with a position and a index
'''

class Node:
    id_iter = itertools.count(start=-2, step=1) #not counting properly
    def __init__(self):
        self.position = np.array([0,0,0])
        self.index = next(self.id_iter)
        self.connectors = []
        self.orientation_node = False
        self.dsDNA = True
    
    def set_position(self, e0):
        self.position = np.asarray(e0)        
    
    def get_position(self):
        return self.position

class OrientationNode(Node):
    def __init__(self, *args):
        self.connectors = []
        self.position = np.array([0, 0, 0])
        self.orientation_vector = np.array([0, 0, 0])
        self.orientation_node = True

class Connector:
    def __init__(self, node_1, node_2):
        assert isinstance(node_1, Node)
        assert isinstance(node_2, Node)
        self.node_1 = node_1
        self.node_2 = node_2
        for node in self.node_1, self.node_2:
            if self not in node.connectors:
                node.connectors.append(self)
                
    @property
    def node_index(self) -> Tuple[int]:
        """Indices of nodes, with node_1 first and node_2 second
        """
        return (self.node_1.index, self.node_2.index)
    @property
    def vector(self) -> np.ndarray:
        """Absolute vector from node_1 to node_2
        """
        return self.node_2.position - self.node_1.position
    @property
    def length(self) -> float:
        """Distance between a connector's two nodes
        """
        return np.linalg.norm(self.vector)
    @property
    def direction(self) -> np.ndarray:
        """Unit vector of direction from node_1 to node_2
        """
        return self.vector / self.length
     
class Beam(Connector):
    def __init__(self, *args):
        super().__init__(*args)
        self.orientation_beam = False
        def find_angle(v1, v2):
            v1_u = v1 / np.linalg.norm(v1)
            v2_u = v2 / np.linalg.norm(v2)
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

if __name__ == "__main__":
    nodes = [
        Node(),
        Node(),
    ]
    connector = Connector(nodes[0], nodes[1])
    for i, node in enumerate(nodes):
        #print(f"Node {i} Connectors: {node.connectors}")
        assert connector in node.connectors   
                    
def set_orientation_nodes(system):
    orientation_elements = []
    for node in system:
        if isinstance(node, Node):
            if node.dsDNA == True:
                temp_node = OrientationNode()
                temp_node.set_position(np.array([node.position[0] + 0.15, node.position[1], node.position[2]])) #0.15 in nm
                temp_beam = Beam(temp_node, node)
                temp_beam.orientation_beam = True
                temp_beam.orientation_vector = node.position - temp_node.position
                orientation_elements.append(temp_node)
                orientation_elements.append(temp_beam)
    system = system + orientation_elements
    return system


def calculate_total_energy(system, T = 298): #T = temperature (K)
    K_B = 1.38065E-23 
    k_bT = K_B * T
    total_energy = 0
    
    def calculate_bond_energy(r1, r2, k, L_0): #pos vectors r1 and r2 of nodes 1 and 2 
        x1, x2, y1, y2, z1, z2 = r1[0], r2[0], r1[1], r2[1], r1[2], r2[2]
        bond_energy = (1/2) * k * ((np.sqrt(np.abs(x1-x2)+np.abs(y1-y2)+np.abs(z1-z2)))-L_0)**2
        return bond_energy
    
    def calculate_angular_energy(r1, r2, r3, k, theta_0):
        v1_u = (r2-r1) / np.linalg.norm(r2-r1)
        v2_u = (r3-r2) / np.linalg.norm(r3-r2)
        theta = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        angular_energy = (1/2) * k * (theta - theta_0) ** 2
        return angular_energy
        
    def calculate_dihedral_energy(r1, r2, r3, r4, k, phi_0):
        pass
    
    # bond energy calculation
    for element in system:
        if isinstance(element, Beam):
            if (element.node_1.orientation_node == False and element.node_2.orientation_node == False) and (element.node_1.dsDNA and element.node_2.dsDNA): 
                # both n1 and n2 are backbone nodes and dsDNA
                L_0 = 0.34 # nm
                k_dsDNA = 2.94 # J/m^2
                energy = calculate_bond_energy(element.node_1.position, element.node_2.position, k_dsDNA, L_0)
                total_energy += energy
            elif element.node_1.orientation_node and element.node_2.dsDNA:
                # Energy of orientation node connected to a real node
                assert (element.node_1.orientation_node and element.node_2.dsDNA), "Only dsDNA have orientation beads"
                L_0 = 0.15 #nm
                k_orient = 20.8 #J/m^2
                energy = calculate_bond_energy(element.node_1.position, element.node_2.position, k_orient, L_0)
                total_energy += energy
            elif (element.node_1.orientation_node == False and element.node_2.orientation_node == False) and (element.node_1.dsDNA == False and element.node_2.dsDNA == False):
                # both nodes are ssDNA 
                L_0 = 0.64 # nm
                k_ssDNA = 2.352 # J/m^2
                energy = calculate_bond_energy(element.node_1.position, element.node_2.position, k_ssDNA, L_0)
                total_energy += energy
            elif (element.node_1.orientation_node == False and element.node_2.orientation_node == False) and (element.node_1.dsDNA == True and element.node_2.dsDNA == False):
                # first node is dsDNA and second node is ssDNA, L_0 and k are averages of the two L_0 and k values for ss and dsDNA
                L_0 = 0.49 # nm
                k_average = (2.94 + 2.352)/2
                energy = calculate_bond_energy(element.node_1.position, element.node_2.position, k_average, L_0)
                total_energy += energy
                
        elif isinstance(element, Node):
            pass

    return total_energy

        

def main():
    return

if __name__ == '__main__':
    pass
    #main()
    
#CanDo system 1 (5 beads connected)
node_1, node_2, node_3, node_4, node_5 = Node(), Node(), Node(), Node(), Node()

node_1.set_position([0, 0, 0])
node_2.set_position([1, 0, 0])
node_3.set_position([2, 0, 0])
node_4.set_position([3, 0, 0])
node_5.set_position([4, 0, 0])

beam_1, beam_2, beam_3, beam_4, = Beam(node_1, node_2), Beam(node_2, node_3), Beam(node_3, node_4), Beam(node_4, node_5)

system_1 = [node_1, node_2, node_3, node_4, node_5, beam_1, beam_2, beam_3, beam_4]
system_1 = set_orientation_nodes(system_1)
print("Energy of system 1 is", calculate_total_energy(system_1))

#System 2
node_6, node_7 = Node(), Node()
node_6.set_position([5, 0, 0])
node_7.set_position([6, 0, 0])
beam_5, beam_6 = Beam(node_5, node_6), Beam(node_6, node_7)
node_4.dsDNA = False
system_2 = [node_1, node_2, node_3, node_4, node_5, node_6, node_7, beam_1, beam_2, beam_3, beam_4, beam_5, beam_6]
system_2 = set_orientation_nodes(system_2)
print("Energy of system 2 is", calculate_total_energy(system_2))
