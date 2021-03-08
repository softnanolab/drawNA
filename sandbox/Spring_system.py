import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import itertools

import pandas as pd
'''
The Node class creates a node object that has a position accessed by self.get_position. The set_position method allows the user
to input a numpy array with elements 0,1,2,3 representing the e0, e1, e2, and e3 vectors, respectively. e0 is the position vector of the node, e1
points towards the major groove, e2 points towards the base in a chosen strand and e3 is the 
vector product of e1 and e2
'''

class Node:
    id_iter = itertools.count()
    def __init__(self):
        self.position = np.array([[0,0,0], [0,0,0], [0,0,0], [0,0,0]])
        self.index = next(Node.id_iter)
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

# test_node = Node()
# test_node.set_position([1,1,1], [0, -1, -1], [0, -2, -1])
# print(type(test_node.position)) # returns ndarray
# print(test_node.position)

class Connector:
    def __init__(self, node_1, node_2):
        assert isinstance(node_1, Node)
        assert isinstance(node_2, Node)
        self.nodes = node_1, node_2
        self.node_index = [node_1.index, node_2.index]
        #Calculates unit vectors of node 1 for e1, e2 and e3 vectors
        self.N1_e1_u = node_1.position[1] / np.linalg.norm(node_1.position[1])
        self.N1_e2_u = node_1.position[2] / np.linalg.norm(node_1.position[2])
        self.N1_e3_u = node_1.position[3] / np.linalg.norm(node_1.position[3])
        #Calculates unit vectors of node 2 for e1, e2 and e3 vectors
        self.N2_e1_u = node_2.position[1] / np.linalg.norm(node_2.position[1])
        self.N2_e2_u = node_2.position[2] / np.linalg.norm(node_2.position[2])
        self.N2_e3_u = node_2.position[3] / np.linalg.norm(node_2.position[3]) 

    def connect_nodes(self, node_1, node_2):
        self.starting_pos = node_1.position[0]
        self.ending_pos = node_2.position[0]
        return self
    def find_angle(self, node_1, node_2):
        def unit_vector(vector):
            return vector / np.linalg.norm(vector)
        v1_u = unit_vector(node_1.position[2])
        v2_u = unit_vector(node_2.position[2])
        self.twist_angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        v3_u = unit_vector(node_1.position[1])
        v4_u = unit_vector(node_2.position[1])
        self.bend_angle_1 = np.arccos(np.clip(np.dot(v3_u, v4_u), -1.0, 1.0))
        v5_u = unit_vector(node_1.position[3])
        v6_u = unit_vector(node_2.position[3])
        self.bend_angle_2 = np.arccos(np.clip(np.dot(v5_u, v6_u), -1.0, 1.0))
        return self

class Beam(Connector):
    def __init__(self, *args):
        super().__init__(*args)
        super().connect_nodes(*args)
        super().find_angle(*args)        
        self.duplex = True
        self.stretch_mod = 1100 #pN
        self.bend_mod = 230 #pN nm^2
        self.twist_mod = 460 #pN nm^2
        self.orientation = args[1].position[0] - args[0].position[0]
        self.length = np.linalg.norm(self.orientation)  
        
    def get_orientation(self):
        return self.orientation
    def get_length(self):
        return self.length
    def get_angles(self):
        return self.bend_angle_1, self.bend_angle_2, self.twist_angle
    def nick_beam(self):
        self.duplex = False
        self.bend_mod = self.bend_mod/100
        self.twist_mod = self.twist_mod/100
        
class Spring(Connector):
    def __init__(self, *args, bulge=False):
        super().__init__(*args)
        self.angle = super().find_angle(*args)
        
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
            beam_start_pos.append(element.starting_pos)
            beam_end_pos.append(element.ending_pos)
                                
    ax = plt.axes(projection = "3d")
    ax.scatter3D(x_pos_node, y_pos_node, z_pos_node, color = "green")
    
    plt.show()

def energy_function(system):
    nodes = []
    beams = []
    springs = []
    for element in system:
        if isinstance(element, Node):
            nodes.append(element)
        elif isinstance(element, Beam):
            beams.append(element)
        elif isinstance(element, Spring):
            springs.append(element)
    print(nodes, "\n", beams, "\n", springs)
    
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
test_beam_1 = Beam(test_node_1, test_node_2)
test_beam_2 = Beam(test_node_2, test_node_3)
test_beam_3 = Beam(test_node_3, test_node_4)
test_beam_4 = Beam(test_node_4, test_node_1)
square_system = [test_node_1, test_node_2, test_node_3, test_beam_4,
          test_beam_1, test_beam_2, test_beam_3, test_node_4
          ]
test_spring = Spring(test_node_1, test_node_2)
test_system = [test_node_1, test_node_2, test_beam_1] 

energy_function(square_system)  
def main():
    matrix = generate_matrix(square_system)
    #visualise_spring_system(system)
    visualise(matrix, show=True)
    return

if __name__ == '__main__':
    pass
    #main()