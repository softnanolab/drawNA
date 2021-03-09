import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import itertools
from typing import Tuple
import pandas as pd
'''
The Node class creates a node object that has a position accessed by self.get_position. The set_position method allows the user
to input a numpy array with elements 0,1,2,3 representing the e0, e1, e2, and e3 vectors, respectively. e0 is the position vector of the node, e1
points towards the major groove, e2 points towards the base in a chosen strand and e3 is the 
vector product of e1 and e2
'''

class Node:
    id_iter = itertools.count(start=0, step=1)
    def __init__(self):
        self.position = np.array([[0,0,0], [0,0,0], [0,0,0], [0,0,0]])
        self.index = next(self.id_iter)-2
        self.connectors = []
    def set_position(self, e0, e1, e2, e3):
        '''
        Parameters: 
        Input: e0, e1, e2, e3 (ndarray)  
        Output: [e0, e1, e2, e3] (ndarray) dim=4x3 
        '''
        self.position[0] = np.asarray(e0)
        self.position[1] = np.asarray(e1)
        self.position[2] = np.asarray(e2)
        self.position[3] = np.asarray(e3)
        
    def get_position(self):
        return self.position

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
        return self.node_2.position[0] - self.node_1.position[0]
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
        self.duplex = True
        self.stretch_mod = 1100 #pN
        self.bend_mod = 230 #pN nm^2
        self.twist_mod = 460 #pN nm^2
        def find_angle(v1, v2):
            v1_u = v1 / np.linalg.norm(v1)
            v2_u = v2 / np.linalg.norm(v2)
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        self.twist_angle = find_angle(self.node_1.position[2], self.node_2.position[2])


    def nick_beam(self):
        self.duplex = False
        self.bend_mod = self.bend_mod/100
        self.twist_mod = self.twist_mod/100

if __name__ == "__main__":
    nodes = [
        Node(),
        Node(),
    ]
    connector = Connector(nodes[0], nodes[1])
    for i, node in enumerate(nodes):
        print(f"Node {i} Connectors: {node.connectors}")
        assert connector in node.connectors   
        
class Spring(Connector):
    def __init__(self, *args, bulge=False):
        super().__init__(*args)
        
        self.rot_x = 1353 #pN nm rad-1
        self.rot_y = 1353 #pN nm rad-1
        self.rot_z = 135.3 #pN nm rad-1

        self._bulge = False

        self.bulge = bulge

    @property
    def bulge(self):
        return self._bulge

    @bulge.setter
    def bulge(self, value: bool):
        assert isinstance(value, bool)
        self._bulge = value
        if value:
            self.rot_x = 13.53 #pN nm rad-1
            self.rot_y = 0 #pN nm rad-1
            self.rot_z = 13.53 #pN nm rad-1
        else:
            self.rot_x = 1353 #pN nm rad-1
            self.rot_y = 1353 #pN nm rad-1
            self.rot_z = 135.3 #pN nm rad-1
            

def visualise_spring_system(system):    
    pass
    node_pos = []
    x_pos_node = []
    y_pos_node = []
    z_pos_node = []
    beam_start_pos = []
    beam_end_pos = []
    
    for element in system:
        if isinstance(element, Node):
            node_pos.append(element.position[0])
            x_pos_node.append(element.position[0][0])
            y_pos_node.append(element.position[0][1])
            z_pos_node.append(element.position[0][2])
        elif isinstance(element, Beam):
            pass
                                
    ax = plt.axes(projection = "3d")
    ax.scatter3D(x_pos_node, y_pos_node, z_pos_node, color = "green")
    
    plt.show()

def energy_function(system):
    energy = 0
    L0 = 0.34 #nm equilibrium bond distance
    eq_theta_twist = 0.5985 # 34.29 degrees
    def find_angle(v1, v2):
        v1_u = v1 / np.linalg.norm(v1)
        v2_u = v2 / np.linalg.norm(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    for element in system:
        if isinstance(element, Beam):
            stretch_energy_x = (1/2) * element.stretch_mod * (np.absolute(element.node_1.position[0][0] - element.node_2.position[0][0])-L0) ** 2
            stretch_energy_y = (1/2) * element.stretch_mod * (np.absolute(element.node_1.position[0][1] - element.node_2.position[0][1])-L0) ** 2
            stretch_energy_z = (1/2) * element.stretch_mod * (np.absolute(element.node_1.position[0][0] - element.node_2.position[0][0])-L0) ** 2
            twist_energy = (1/2) * element.twist_mod * (element.twist_angle - eq_theta_twist) ** 2
            energy += stretch_energy_x + stretch_energy_y + stretch_energy_z + twist_energy
        elif isinstance(element, Node):
            if len(element.connectors) == 2:
                bend_angle =  find_angle(element.connectors[0].direction, element.connectors[1].direction)
                bend_energy = (1/2) * element.connectors[0].bend_mod * (bend_angle - 0) **2
                energy += bend_energy
                
    return energy
def generate_matrix(system):  
    #print(pd.DataFrame(matrix))        
    return -1
   
def visualise(matrix: np.ndarray, show: bool = False):
    """Taken from here: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
    """
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.imshow(matrix)
    # Loop over data dimensions and create text annotations.
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            ax.text(j, i, matrix[i, j], ha="center", va="center", color="w", fontsize='xx-small')
    # plot lines
    n_elements = len(matrix) // 6
    for i in range(n_elements - 1):
        pos = (i+1) * 6 - 0.5
        ax.plot([-0.5, len(matrix)-0.5], [pos, pos], 'k-')
        ax.plot([pos, pos], [-0.5, len(matrix)-0.5], 'k-')
    if show:
        plt.show()

test_node_1 = Node()
test_node_1.set_position([0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0])

test_node_2 = Node()
test_node_2.set_position([1, 0, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0])

test_node_3 = Node()
test_node_3.set_position([1, 1, 0], [2, 1, 1], [1,-1,-2], [1, 1, 1])

test_node_4 = Node()
test_node_4.set_position([0, 1, 0], [1, 1, 1], [0, 2, 2], [1, 1, 1])

# Need to assert that only the elements used are initialised
test_beam_1 = Beam(test_node_1, test_node_2)
test_beam_2 = Beam(test_node_2, test_node_3)

'''
test_beam_3 = Beam(test_node_3, test_node_4)
test_beam_4 = Beam(test_node_4, test_node_1)



square_system = [test_node_1, test_node_2, test_node_3, test_node_4, 
          test_beam_1, test_beam_2, test_beam_3, test_beam_4,
          ]

'''
test_system = [test_node_1, test_node_2, test_beam_1] 
print("Energy = ", energy_function(test_system))

def main():    
    #visualise_spring_system(system)
    visualise(matrix, show=True)
    return

if __name__ == '__main__':
    pass
    #main()