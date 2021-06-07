import numpy as np
import itertools
import scipy.optimize
from typing import Tuple
import matplotlib.pyplot as plt
from drawNA.readers import OXDNAReader
import random

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
            r1 = self.nodes[i].position
            r2 = self.nodes[i+1].position
            r3 = self.nodes[i+2].position
            
            
            v1 = r2 - r1
            v2 = r3 - r2

            v1 = v1[0]
            v2 = v2[0]            

            dot_product = (v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2])
            mod = np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2) * np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
            
            x = dot_product/mod
            x = round(x, 4)
            theta = np.arccos(x)
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
                r2 = self.nodes[i].position
                r3 = self.nodes[i+1].position
                r2 = r2[0]
                r3 = r3[0]
                
                r1 = np.array([r2[0]+0.15, r2[1], r2[2]])
                r4 = np.array([r3[0]+0.15, r3[1], r3[2]])
                
                r_12 = r2 - r1
                r_23 = r3 - r2
                r_34 = r4 - r3
                
                v1_dihedral = np.cross(r_12, r_23)
                v2_dihedral = np.cross(r_23, r_34)
                
                v1_mod = np.sqrt((v1_dihedral[0]**2) + (v1_dihedral[1]**2) + (v1_dihedral[2]**2))
                v2_mod = np.sqrt((v2_dihedral[0]**2) + (v2_dihedral[1]**2) + (v2_dihedral[2]**2))
                
                v1_dihedral_u = v1_dihedral/v1_mod
                v2_dihedral_u = v2_dihedral/v2_mod
                
                dot_product = (v1_dihedral_u[0]*v2_dihedral_u[0]) + (v1_dihedral_u[1]*v2_dihedral_u[1]) + (v1_dihedral_u[2]*v2_dihedral_u[2])
                dot_product = round(dot_product, 7)
                
                phi = np.arccos(dot_product)
                phi_0 = 0.602
                k = 132.35 * kbT                
                total_energy += 0.5 * k * (phi - phi_0)**2
            # Improper dihedral
            if self.nodes[i].dsDNA and self.nodes[i+1].dsDNA and self.nodes[i+2].dsDNA:
                r1 = self.nodes[i].position
                r2 = self.nodes[i+1].position
                r1 = r1[0]
                r2 = r2[0]
                r3 = np.array([r2[0]+0.15, r2[1], r2[2]])
                r4 = self.nodes[i+2].position

                r_12 = r2-r1
                r_23 = r3-r2
                r_32 = r2-r3
                r_24 = r4-r2
                
                
                v1_dihedral = np.cross(r_12, r_23)
                v2_dihedral = np.cross(r_32, r_24)
                v2_dihedral = v2_dihedral[0]
                
                v1_mod = np.sqrt((v1_dihedral[0]**2) + (v1_dihedral[1]**2) + (v1_dihedral[2]**2))
                v2_mod = np.sqrt((v2_dihedral[0]**2) + (v2_dihedral[1]**2) + (v2_dihedral[2]**2))
                
                v1_dihedral_u = (v1_dihedral)/v1_mod
                v2_dihedral_u = (v2_dihedral)/v2_mod
                
                dot_product = (v1_dihedral_u[0]*v2_dihedral_u[0]) + (v1_dihedral_u[1]*v2_dihedral_u[1]) + (v1_dihedral_u[2]*v2_dihedral_u[2])
                
                phi = np.arccos(dot_product)
                phi_0 = 0
                k = 13.755 * kbT
                
                total_energy += 0.5 * k * (phi-phi_0)**2

        # Orientation node contributions
        # Stretching contributions:
        for i in range(1, len(self.nodes)):
            if self.nodes[i].dsDNA:
                r1 = self.nodes[i].position
                r1 = r1[0]
                r2 = np.array([r1[0]+0.15, r1[1], r1[2]])
                k = 20.8
                L_0 = 0.15
                total_energy += 0.5 * k * (np.linalg.norm(r2-r1) - L_0)**2
        # Bending contributions
        for i in range(len(self.nodes)-1):
            if self.nodes[i].dsDNA:
                r2 = self.nodes[i].position
                r2 = r2[0]
                r1 = np.array([r2[0]+0.15, r2[1], r2[2]])
                r3 = self.nodes[i+1].position
                r3 = r3[0]
                k = 13.755 * kbT
                theta_0 = np.pi/2
          
                v1 = np.asarray((r2 - r1))
                v2 = np.asarray((r3 - r2))
                dot_product = (v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2])
                mod = np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2) * np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                
                x = dot_product / mod
                x = round(x, 8)
                theta = np.arccos(x)
                
                total_energy += 0.5 * k * (theta-theta_0)**2
        
        # Crossovers (type B)
        for i in range(len(self.connector)):
            if self.connector[i].node_index[1] != self.connector[i].node_index[0]+1:
                r1 = self.nodes[self.connector[i].node_index[0]].position
                r2 = self.nodes[self.connector[i].node_index[1]].position
                r3 = self.nodes[self.connector[i].node_index[1]-1].position
                r1 = r1[0]
                r2 = r2[0]
                r3 = r3[0]
                
                r4 = np.array([r2[0]+0.15, r2[1], r2[2]])
                r_12 = r2-r1
                r_23 = r3-r2
                r_24 = r4-r2
                
                v1_crossover = np.cross(r_12, r_23)
                v2_crossover = np.cross(r_23, r_24)
                
                v1_mod = np.sqrt((v1_crossover[0]**2) + (v1_crossover[1]**2) + (v1_crossover[2]**2))
                v2_mod = np.sqrt((v2_crossover[0]**2) + (v2_crossover[1]**2) + (v2_crossover[2]**2))
                
                v1_crossover_u = v1_crossover/v1_mod
                v2_crossover_u = v2_crossover/v2_mod
                
                dot_product = (v1_crossover_u[0]*v2_crossover_u[0]) + (v1_crossover_u[1]*v2_crossover_u[1]) + (v1_crossover_u[2]*v2_crossover_u[2])
                dot_product = round(dot_product, 8)
                phi = np.arccos(dot_product)
                phi_0 = (2*np.pi)/3
                k = 13.755 * kbT
                
                total_energy += 0.5 * k * (phi-phi_0)**2
            # Crossovers (type A)
            if self.connector[i].node_index[1] != self.connector[i].node_index[0]+1:
                r1 = self.nodes[self.connector[i].node_index[0]-1].position
                r2 = self.nodes[self.connector[i].node_index[0]].position
                r3 = self.nodes[self.connector[i].node_index[1]].position
                r4 = self.nodes[self.connector[i].node_index[1]+1].position
                
                r1 = r1[0]
                r2 = r2[0]
                r3 = r3[0]
                r4 = r4[0]
                
                r_12 = r2-r1
                r_23 = r3-r2
                r_34 = r4-r3
                
                v1_crossover = np.cross(r_12, r_23)
                v2_crossover = np.cross(r_23, r_24)
                
                v1_mod = np.sqrt((v1_crossover[0]**2) + (v1_crossover[1]**2) + (v1_crossover[2]**2))
                v2_mod = np.sqrt((v2_crossover[0]**2) + (v2_crossover[1]**2) + (v2_crossover[2]**2))
                v1_crossover_u = v1_crossover/v1_mod
                v2_crossover_u = v2_crossover/v2_mod
                
                
                dot_product = (v1_crossover_u[0]*v2_crossover_u[0]) + (v1_crossover_u[1]*v2_crossover_u[1]) + (v1_crossover_u[2]*v2_crossover_u[2])
                dot_product = round(dot_product, 8)
                phi = np.arccos(dot_product)
                
                k = 13.755 * kbT
                phi_0 = 0
                total_energy += 0.5 * k * (phi-phi_0)**2
            
            #stretching crossover
            if self.connector[i].node_index[1] != self.connector[i].node_index[0]+1:
                r1 = self.nodes[self.connector[i].node_index[0]].position[0]
                r2 = self.nodes[self.connector[i].node_index[1]].position[0]
                l_0 = 1.85
                k = 2.77 
                l = np.linalg.norm(r2-r1)
                total_energy += 0.5 * k * (l-l_0)**2
        
        
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
        coordinates = []
        for i in range((len(solution.x)//3)):
            node_pos = []
            node_pos.append(solution.x[3*i])
            node_pos.append(solution.x[3*i+1])
            node_pos.append(solution.x[3*i+2])
            coordinates.append(node_pos)
        return coordinates
    
    def visualise(self):
        nodes = self.nodes
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xs = []
        ys = []
        zs = []
        for i in range(len(nodes)):
            xs.append(nodes[i].position[0])
            ys.append(nodes[i].position[1])
            zs.append(nodes[i].position[2])
        for i in range(len(xs)-1):    
            ax.quiver(xs[i], ys[i], zs[i], xs[i+1]-xs[i], ys[i+1]-ys[i], zs[i+1]-zs[i])
            pass
        for i in range(len(self.connector)):
            
            if self.connector[i].node_index[1] != self.connector[i].node_index[0]+1:
                node_1_pos = self.nodes[self.connector[i].node_index[0]].position
                node_2_pos = self.nodes[self.connector[i].node_index[1]].position
                     
                
                ax.quiver(node_1_pos[0], node_1_pos[1], node_1_pos[2],
                          (node_2_pos[0]-node_1_pos[0]), (node_2_pos[1]-node_1_pos[1]),
                          (node_2_pos[2]-node_1_pos[2]))
                
        
        ax.scatter(xs, ys, zs, label="Non-equilibrium structure")


        node_positions = self.minimise()       
        xs = []
        ys = []
        zs = []
        for i in range(len(node_positions)):
            xs.append(node_positions[i][0])
            ys.append(node_positions[i][1])
            zs.append(node_positions[i][2])
        
        for i in range(len(xs)-1):
            #print("Vectors are: ", [xs[i+1]-xs[i], ys[i+1]-ys[i], zs[i+1]-zs[i]])
            #print("Vectors are: ", [xs[i], ys[i], zs[i]])
            pass
        ax.scatter(xs, ys, zs, label="Equilibrium structure", color="tab:orange")
        
        for i in range(len(xs)-1):
            ax.quiver(xs[i], ys[i], zs[i], xs[i+1]-xs[i], ys[i+1]-ys[i], zs[i+1]-zs[i], color="tab:orange", length=1)
            pass
        for i in range(len(self.connector)):
            if self.connector[i].node_index[1] != self.connector[i].node_index[0]+1:
                node_1_pos = list(self.nodes[self.connector[i].node_index[0]].position[0])
                node_2_pos = list(self.nodes[self.connector[i].node_index[1]].position[0])
                
                
                ax.quiver(node_1_pos[0], node_1_pos[1], node_1_pos[2],
                          node_2_pos[0]-node_1_pos[0], node_2_pos[1]-node_1_pos[1],
                          node_2_pos[2]-node_1_pos[2], color="tab:orange")  
                
                
                pass
        
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_zticklabels([])
                
        #Set the scale of the graph
        ax.set_xlim3d(-0.5, 4.5)
        ax.set_ylim3d(-0.5, 4.5)
        ax.set_zlim3d(-0.5, 0.5)
        
        ax.legend()
        
        plt.show()

def reader(file_path):
    f = file_path
    system = OXDNAReader(f).system
    nucleotides = system.nucleotides

    nodes = []
    connectors = []
    
    for i in range(len(nucleotides)):
        node_pos = nucleotides[i].pos_com
        #workaround for a bug in the minimiser
        
        random_displacement = np.array([random.randint(0,9)*0.02, \
        random.randint(0,9)*0.02, random.randint(0,9)*0.02]) 
        node_pos_fixed = node_pos + random_displacement 
        
        
        nodes.append(Node(node_pos))
        
    for i in range(len(nodes)-1):
        connectors.append(Bond(nodes[i], nodes[i+1]))
    
    system = System(nodes, connectors)
    
    return system

def main():
    return

def generate_lattice(dim):
      node_pos = []
      starting_pos = []
      random.seed(42)      
      for i in range(dim):
        array = np.array([0., 0.5*i, 0.])
        starting_pos.append(array)
      for j in range(dim):
          for k in range(dim):
              if j % 2 == 0:
                  pos = starting_pos[j] + np.array([0.5*k, 0., 0.])
                  random_displacement = np.array([random.randint(1, 50)*0.001, random.randint(1, 50)*0.001, random.randint(1, 50)*0.001])
                  pos = pos + random_displacement
                  node_pos.append(list(pos))
              elif j % 2 != 0:
                  pos = starting_pos[j] + np.array([0.5*dim, 0., 0.]) - np.array([0.5*(k+1), 0., 0.])
                  random_displacement = np.array([random.randint(1, 50)*0.001, random.randint(1, 50)*0.001, random.randint(1, 50)*0.001])
                  pos = pos + random_displacement
                  node_pos.append(list(pos))
     
      return node_pos

def generate_system(node_pos):
    nodes = []
    connectors = []
    for i in range(len(node_pos)):
        node = Node(node_pos[i])
        nodes.append(node)
    for i in range(len(nodes)-1):
        bond = Bond(nodes[i], nodes[i+1])
        connectors.append(bond)
        bond_0 = Bond(nodes[1], nodes[12])
        bond_1 = Bond(nodes[17], nodes[24])
        bond_2 = Bond(nodes[33], nodes[38])
        
        connectors.append(bond_0)
        connectors.append(bond_1)
        connectors.append(bond_2)
        
    system = System(nodes, connectors)
    return system

lattice = generate_lattice(7)
example_1 = generate_system(lattice)
example_1.visualise()


