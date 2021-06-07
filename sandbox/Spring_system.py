import numpy as np
import itertools
import scipy.optimize
from typing import Tuple

class Node():
    id_iter = itertools.count(start=-2, step=1) #not counting properly
    def __init__(self, position=[0,0,0], dsDNA=True):
        self.position = np.asarray(position)
        self.index = next(self.id_iter)
        self.connectors = []
        self.dsDNA = dsDNA
    
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
    
class Bond(Connector):
    def __init__(self, *args):
        super().__init__(*args)

if __name__ == "__main__":
    nodes = [
        Node(),
        Node(),
    ]
    connector = Connector(nodes[0], nodes[1])
    for i, node in enumerate(nodes):
        #print(f"Node {i} Connectors: {node.connectors}")
        assert connector in node.connectors   

class System():
    def __init__(self, nodes, connectors, T=298):
        self.nodes = np.asarray(nodes)
        self.connector = np.asarray(connectors)
        self.T = T

    def to_array(self, flatten=True):
        array = []
        for node in self.nodes:
            array.append(node.position)
        array = np.asarray(array)
        return array
   
    def from_array(self, array):                
        for i in range(len(self.nodes)):
            position = []
            position.append([array[3*i], array[3*i+1], array[3*i+2]])
            position = np.asarray(position)
            self.nodes[i].position = position            
        return
    
    @property
    def energy(self):
        total_energy = 0
        T = self.T
        kbT = 1.38065E-5 * T
        
        # Stretching energy
        for i in range(len(self.nodes)-1):
            r1 = self.nodes[i].position[0]
            r2 = self.nodes[i+1].position[0]
            if self.nodes[i].dsDNA and self.nodes[i+1].dsDNA:
                k = 2.94
                L_0 = 0.34 #nm
                total_energy += 0.5 * k * (np.linalg.norm(r2-r1) - L_0)**2
            elif self.nodes[i].dsDNA == False and self.nodes[i+1].dsDNA == False:
                k = 2.352
                L_0 = 0.64
                total_energy += 0.5 * k * (np.linalg.norm(r2-r1) - L_0)**2
            elif (self.nodes[i].dsDNA and self.nodes[i+1].dsDNA == False) or (self.nodes[i].dsDNA == False and self.nodes[i+1].dsDNA):
                k = np.sqrt(2.94 * 2.352) # geometric mean = 2.63
                L_0 = (0.34 + 0.64) / 2 # L_0 mean = 0.49
                total_energy += 0.5 * k * (np.linalg.norm(r2-r1) - L_0)**2
            
        # Bending energy
        for i in range(len(self.nodes)-2):
            r1 = self.nodes[i].position[0]
            r2 = self.nodes[i+1].position[0]
            r3 = self.nodes[i+2].position[0]          
            v1 = np.asarray((r2 - r1))
            v2 = np.asarray((r3 - r2))
            dot_product = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
            mod = np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2) + np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
            theta = np.arccos(dot_product / mod)

            if self.nodes[i].dsDNA and self.nodes[i+1].dsDNA and self.nodes[i+2].dsDNA:
                k = 110.2 * kbT
                theta_0 = 0
                total_energy += 0.5 * k * (theta-theta_0)**2

            elif self.nodes[i].dsDNA == False and self.nodes[i+1] == False and self.nodes[i+2] == False:
                k = 1.675 * kbT
                theta_0 = 0
                total_energy += 0.5 * k * (theta - theta_0)**2
            
            elif (self.nodes[i].dsDNA and self.nodes[i+1].dsDNA == False) or (self.nodes[i].dsDNA == False and self.nodes[i+1].dsDNA):
                k = 1.675 * kbT
                theta_0 = 0
                total_energy += 0.5 * k * (theta - theta_0)**2
            
        # Twisting energy
            if self.nodes[i].dsDNA and self.nodes[i+1].dsDNA:
                r2 = nodes[i].position[0]
                r3 = nodes[i+1].position[0]
                r1 = np.array([r2[0]+0.15, r2[1], r2[2]])
                r4 = np.array([r3[0]+0.15, r3[1], r3[2]])
                
                r_12 = r2 - r1
                r_23 = r3 - r2
                r_34 = r4 - r3
                
                v1_dihedral = np.cross(r_12, r_23)
                v2_dihedral = np.cross(r_23, r_34)

                dot_product = v1_dihedral[0]*v2_dihedral[0]+v1_dihedral[1]*v2_dihedral[1]+v1_dihedral[2]*v2_dihedral[2]
                mod = np.sqrt(v1_dihedral[0]**2 + v1_dihedral[1]**2 + v1_dihedral[2]**2) + np.sqrt(v2_dihedral[0]**2 + v2_dihedral[1]**2 + v2_dihedral[2]**2)
                phi = np.arccos(dot_product / mod)
                phi_0 = 0.602
                k = 132.35 * kbT

                total_energy += 0.5 * k * (phi - phi_0)**2
            # Improper dihedral
            if self.nodes[i].dsDNA and self.nodes[i+1].dsDNA and self.nodes[i+2].dsDNA:
                r1 = self.nodes[i].position[0]
                r2 = self.nodes[i+1].position[0]
                r3 = np.array([r2[0]+0.15, r2[1], r2[2]])
                r4 = self.nodes[i+2].position[0]

                r_12 = r2-r1
                r_23 = r3-r2
                r_32 = r2-r3
                r_24 = r4-r2

                v1_dihedral = np.cross(r_12, r_23)
                v2_dihedral = np.cross(r_32, r_24)

                dot_product = v1_dihedral[0]*v2_dihedral[0]+v1_dihedral[1]*v2_dihedral[1]+v1_dihedral[2]*v2_dihedral[2]
                mod = np.sqrt(v1_dihedral[0]**2 + v1_dihedral[1]**2 + v1_dihedral[2]**2) + np.sqrt(v2_dihedral[0]**2 + v2_dihedral[1]**2 + v2_dihedral[2]**2)
                phi = np.arccos(dot_product / mod)
                phi_0 = 0
                k = 13.755 * kbT
                
                total_energy += 0.5 * k * (phi-phi_0)**2

        # Orientation node contributions
        # Stretching contributions:
        for i in range(len(self.nodes)):
            r1 = self.nodes[i].position[0]
            r2 = np.array([r1[0]+0.15, r1[1], r1[2]])
            k = 20.8
            L_0 = 0.15
            total_energy += 0.5 * k * (np.linalg.norm(r2-r1) - L_0)**2
        # Bending contributions
        for i in range(len(self.nodes)-1):
            if self.nodes[i].dsDNA:
                r2 = self.nodes[i].position[0]
                r1 = np.array([r2[0]+0.15, r2[1], r2[2]])
                r3 = self.nodes[i+1].position[0]
                k = 13.755 * kbT
                theta_0 = np.pi/2
          
                v1 = np.asarray((r2 - r1))
                v2 = np.asarray((r3 - r2))
                dot_product = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
                mod = np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2) + np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                theta = np.arccos(dot_product / mod)
                total_energy += 0.5 * k * (theta-theta_0)**2

        return total_energy
    
    def minimise(self, **kwargs):
        def _e(pos):
            self.from_array(pos)
            return self.energy
        initial_pos = self.to_array()
        solution = scipy.optimize.minimize(_e, initial_pos, **kwargs)
       
        bond_lengths = []
        for i in range((len(solution.x)//4)):
            node_1 = np.array([solution.x[3*i], solution.x[3*i+1], solution.x[3*i+2]])
            node_2 = np.array([solution.x[3*i+3], solution.x[3*i+4], solution.x[3*i+5]])
            bond_lengths.append(np.linalg.norm(node_2-node_1))
        print(bond_lengths)  
        
        return 
    
    
node_0, node_1, node_2, node_3 = Node([0.2, 0.25, 0.05]), Node([0.4, 0.5, 0.2]), Node([0.7, 0.9, 0.5]), Node([1.5, 1.5, 1.])
nodes = [node_0, node_1, node_2, node_3]
bond_0, bond_1, bond_2 = Bond(node_0, node_1), Bond(node_1, node_2), Bond(node_2, node_3)
connectors = [bond_0, bond_1, bond_2]

example_1 = System(nodes, connectors)
example_1.minimise()


node_0.dsDNA, node_1.dsDNA, node_2.dsDNA, node_3.dsDNA = False, False, False, False

ssDNA_example_1 = System(nodes, connectors) 
ssDNA_example_1.minimise()

node_1.dsDNA, node_2.dsDNA = True, True

mixed_example_1 = System(nodes, connectors)
mixed_example_1.minimise()


def main():
    return


