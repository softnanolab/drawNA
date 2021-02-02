import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import itertools
'''
The Node class creates a node object that has a position accessed by self.get_position. The set_position method allows the user
to input a numpy array with elements 0,1,2 representing the e0, e1 and e2 vectors. e0 is the position vector of the node, e1
points towards the major groove, e2 points towards the base in a chosen strand and e3 is automatically calculated by taking the 
vector product of e1 and e2
'''

class Node:
    id_iter = itertools.count()
    def __init__(self):
        self.position = np.array([[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]])
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
        return self.positon

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
        
    def get_vector(self):
        return self.vector
    def get_length(self):
        return self.length
    def get_angles(self):
        return self.bend_angle_1, self.bend_angle_2, self.twist_angle
    def calculate_energy(self):
        energy_x = 0.5 * self.stretch_mod * (self.ending_pos[0] - self.starting_pos[0])**2
        energy_y = 0.5 * self.stretch_mod * (self.ending_pos[1] - self.starting_pos[1])**2
        energy_z = 0.5 * self.stretch_mod * (self.ending_pos[2] - self.starting_pos[2])**2 
        bend_energy_1 = 0.5 * self.bend_mod * self.bend_angle_1 **2
        bend_energy_2 = 0.5 * self.bend_mod * self.bend_angle_2 ** 2
        twist_energy = 0.5 * self.twist_mod * self.twist_angle ** 2
        total_energy = energy_x + energy_y + energy_z + bend_energy_1 + bend_energy_2 + twist_energy
        return total_energy
    def nick_beam(self):
        self.duplex = False
        self.bend_mod = self.bend_mod/100
        self.twist_mod = self.twist_mod/100
        

class Spring(Connector):
    def __init__(self, *args):
        super().__init__(*args)
        self.extension = True
    def torsional_spring(self, *args):
        self.extension = False
        self.angle = super().find_angle(*args)
        
    
            
test_node_1 = Node()
test_node_1.set_position([0, 0, 0], [0, -1, -1], [0, -2, -1], [1, 1, 1])
test_node_2 = Node()
test_node_2.set_position([0, 1, 0], [1, 1, 1], [2, 2, 1], [1, 1, 1])
test_node_3 = Node()
test_node_3.set_position([1, 1, 0], [2, 1, 1], [1,-1,-2], [1, 1, 1])
test_node_4 = Node()
test_node_4.set_position([1, 0, 0], [1, 1, 1], [0, 2, 2], [1, 1, 1])
test_beam_1 = Beam(test_node_1, test_node_2)
test_beam_2 = Beam(test_node_2, test_node_3)
test_beam_3 = Beam(test_node_3, test_node_4)
test_beam_4 = Beam(test_node_4, test_node_1)
system = [test_node_1, test_node_2, test_node_3, test_beam_4,
          test_beam_1, test_beam_2, test_beam_3, test_node_4
          ]
test_spring = Spring(test_node_1, test_node_2)
test_system = [test_node_1, test_node_2, test_beam_1]

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

def generate_matrix(system):
    #initialize an array of 6n x 6n zeros (make sure there's no double count)
    matrix = np.zeros((len(system)*6,len(system)*6))
    #if element is a beam then there will be stretch+bend+twist, cross terms -2 moduli
    
    #if element is a torsional spring then only twist
    #if element is a nick there will be reduced stretch, bend and twist moduli
    #single or double crossovers, open nicks and bulges modeled as torsional springs
    for element in system:
        pass
generate_matrix(system)
    
        




