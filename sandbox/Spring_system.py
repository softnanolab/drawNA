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
            
    def calculate_energy(self):
        rot_x_energy = 0.5 * self.rot_x * self.bend_angle_1 ** 2
        rot_y_energy = 0.5 * self.rot_y * self.bend_angle_2 ** 2
        rot_z_energy = 0.5 * self.rot_z * self.twist_angle ** 2
        total_energy = rot_x_energy + rot_y_energy + rot_z_energy
        return total_energy
 
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
    node_counter = 0 
    #initialize an array of 6n x 6n zeros (make sure there's no double count)
    for element in system:
        if isinstance(element, Node):
            node_counter += 1
    matrix = np.zeros((12*node_counter,12*node_counter))
    #if element is a beam then there will be stretch+bend+twist, cross terms -2 moduli
    for element in system:
        # generate matrix
        if isinstance(element, Beam):
        # initialise variables
            # e1 vector node 1
            e1_x1 = element.N1_e1_u[0] 
            e1_y1 = element.N1_e1_u[1]
            e1_z1 = element.N1_e1_u[2]
            # e1 vector node 2    
            e1_x2 = element.N2_e1_u[0] 
            e1_y2 = element.N2_e1_u[1]
            e1_z2 = element.N2_e1_u[2]
            # e2 vector node 1 
            e2_x1 = element.N1_e2_u[0] 
            e2_y1 = element.N1_e2_u[1]
            e2_z1 = element.N1_e2_u[2]
            # e2 vector node 2
            e2_x2 = element.N2_e2_u[0] 
            e2_y2 = element.N2_e2_u[1]
            e2_z2 = element.N2_e2_u[2]
            # e3 vector node 1
            e3_x1 = element.N1_e3_u[0] 
            e3_y1 = element.N1_e3_u[1]
            e3_z1 = element.N1_e3_u[2]
            # e3 vector node 2
            e3_x2 = element.N2_e3_u[0] 
            e3_y2 = element.N2_e3_u[1]
            e3_z2 = element.N2_e3_u[2]
            #Useful terms
            e1e1_sum = e1_x1*e1_x2 + e1_y1*e1_y2 + e1_z1*e1_z2
            e2e2_sum = e2_x1*e2_x2 + e2_y1*e2_y2 + e2_z1*e2_z2
            e3e3_sum = e3_x1*e3_x2 + e3_y1*e3_y2 + e3_z1*e3_z2
        
            denom_e1e1_1 = np.sqrt(1-e1e1_sum**2)
            denom_e1e1_2 = 2*(1-(e1e1_sum)**2)**(3/2)
            denom_e1e1_3 = 1-(e1e1_sum)**2
        
            denom_e2e2_1 = np.sqrt(1-e2e2_sum**2)
            denom_e2e2_2 = 2*(1-(e2e2_sum)**2)**(3/2)
            denom_e2e2_3 = 1-(e2e2_sum)**2
        
            denom_e3e3_1 = np.sqrt(1-e3e3_sum**2)
            denom_e3e3_2 = 2*(1-(e3e3_sum)**2)**(3/2)
            denom_e3e3_3 = 1-(e3e3_sum)**2
        
     #Calculate derivatives noting that modulus of unit vectors is 1
        #e1-e1 vector derivatives
            dx1dx1_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(-2*e1_x2*e1e1_sum + 2*e1_x1*(e1e1_sum)**2)*(e1_x2-e1_x1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x2-e1_x1*e1e1_sum)**2)/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (-2*e1_x1*e1_x2 + 3*(e1_x1**2)*e1e1_sum - e1e1_sum)) / denom_e1e1_1)
            
            dx1dx2_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(-2*e1_x1*e1e1_sum + 2*e1_x2*(e1e1_sum)**2)*(e1_x2-e1_x1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x1-e1_x2*e1e1_sum) * (e1_x2 - e1_x1*e1e1_sum))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (e1_x1*e1_x2*e1e1_sum + 1 - e1_x2**2 - e1_x1**2)) / denom_e1e1_1)
            
            dx2dx2_e1 = ((((np.arccos(e1e1_sum))*(element.bend_mod)*(-2*e1_x1*e1e1_sum + 2*e1_x2*(e1e1_sum)**2)*(e1_x1-e1_x2*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x1-e1_x2*e1e1_sum)**2)/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (-2*e1_x1*e1_x2 + 3*(e1_x2**2)*e1e1_sum - e1e1_sum)) / denom_e1e1_1))
            
            dx1dy1_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(-2*e1_y2*e1e1_sum + 2*e1_y1*(e1e1_sum)**2)*(e1_x2-e1_x1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x2-e1_x1*e1e1_sum) * (e1_y2-e1_y1*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (3*e1_x1*e2_y1*e1e1_sum - e1_x1*e1_y2 - e1_x2*e1_y1) / denom_e1e1_1))
            
            dx2dy1_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(-2*e1_y2*e1e1_sum + 2*e1_y1*(e1e1_sum)**2)*(e1_x1-e1_x2*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x1-e1_x2*e1e1_sum) * (e1_y2-e1_y1*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (e1_x2*e1_y1*e1e1_sum - e1_x2*e1_y2 - e1_x1*e1_y1) / denom_e1e1_1))
 
            dx1dy2_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_y2*e1e1_sum**2 + -2*e1_y1*(e1e1_sum))*(e1_x2-e1_x1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x2-e1_x1*e1e1_sum) * (e1_y1-e1_y2*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (e1_x1*e1_y2*e1e1_sum - e1_x2*e1_y2 - e1_x1*e1_y1) / denom_e1e1_1))
            
            dx2dy2_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_y2*e1e1_sum**2 + -2*e1_y1*(e1e1_sum))*(e1_x1-e1_x2*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x1-e1_x2*e1e1_sum) * (e1_y1-e1_y2*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (3*e1_x2*e1_y2*e1e1_sum - e1_x2*e1_y2 - e1_x1*e1_y1) / denom_e1e1_1))
 
            dx1dz1_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_z1*e1e1_sum**2 + -2*e1_z2*(e1e1_sum))*(e1_x2-e1_x1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x2-e1_x1*e1e1_sum) * (e1_z2-e1_z1*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (3*e1_x1*e1_z1*e1e1_sum - e1_x1*e1_z2 - e1_x2*e1_z1) / denom_e1e1_1))
 
            dx2dz1_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_z1*e1e1_sum**2 + -2*e1_z2*(e1e1_sum))*(e1_x1-e1_x2*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_x1-e1_x2*e1e1_sum) * (e1_z2-e1_z1*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (e1_x2*e1_z1*e1e1_sum - e1_x1*e1_z1 - e1_x2*e1_z2) / denom_e1e1_1))
 
            dx1dz2_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_z2*e1e1_sum**2 + -2*e1_z1*(e1e1_sum))*(e1_x2-e1_x1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_z1-e1_z2*e1e1_sum) * (e1_x2-e1_x1*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (e1_x1*e1_z2*e1e1_sum - e1_x1*e1_z1 - e1_x2*e1_z2) / denom_e1e1_1))
 
            dx2dz2_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_z2*e1e1_sum**2 + -2*e1_z1*(e1e1_sum))*(e1_x1-e1_x2*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_z1-e1_z2*e1e1_sum) * (e1_x1-e1_x2*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (3*e1_x2*e1_z2*e1e1_sum - e1_x2*e1_z1 - e1_x1*e1_z2) / denom_e1e1_1))
 
            dy1dy1_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_y1*e1e1_sum**2 + -2*e1_y2*(e1e1_sum))*(e1_y2-e1_y1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_y2-e1_y1*e1e1_sum)**2)/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (-2*e1_y1*e1_y2 + 3*(e1_y1**2)*e1e1_sum - e1e1_sum)) / denom_e1e1_1)
 
            dy1dy2_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_y2*e1e1_sum**2 + -2*e1_y1*(e1e1_sum))*(e1_y2-e1_y1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_y1-e1_y2*e1e1_sum) * (e1_y2-e1_y1*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (e1_y1*e1_y2*e1e1_sum - e1_y1**2 - e1_y2**2 + 1) / denom_e1e1_1))
 
            dy2dy2_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_y2*e1e1_sum**2 + -2*e1_y1*(e1e1_sum))*(e1_y1-e1_y2*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_y1-e1_y2*e1e1_sum)**2)/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (-2*e1_y1*e1_y2 + 3*(e1_y2**2)*e1e1_sum - e1e1_sum)) / denom_e1e1_1)
 
            dy1dz1_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_z1*e1e1_sum**2 + -2*e1_z2*(e1e1_sum))*(e1_y2-e1_y1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_y2-e1_y1*e1e1_sum) * (e1_z2-e1_z1*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (3*e1_y1*e1_z1*e1e1_sum - e1_y2*e1_z1 - e1_y1*e1_z2) / denom_e1e1_1))
 
            dy2dz1_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_z1*e1e1_sum**2 + -2*e1_z2*(e1e1_sum))*(e1_y1-e1_y2*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_y1-e1_y2*e1e1_sum) * (e1_z2-e1_z1*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (e1_y2*e1_z1*e1e1_sum - e1_y2*e1_z2 - e1_y1*e1_z1) / denom_e1e1_1))
 
            dy1dz2_e1 = (((np.arccos(e1e1_sum))*(element.bend_mod)*(2*e1_z2*e1e1_sum**2 + -2*e1_z1*(e1e1_sum))*(e1_y2-e1_y1*e1e1_sum))/denom_e1e1_2
            + ((element.bend_mod) * (e1_y2-e1_y1*e1e1_sum) * (e1_z1-e1_z2*(e1e1_sum)))/denom_e1e1_3
            - ((np.arccos(e1e1_sum)) * (element.bend_mod) * (e1_y1*e1_z2*e1e1_sum - e1_y2*e1_z2 - e1_y1*e1_z1) / denom_e1e1_1))
            
 
     #12 DoF (3 translational 2 bend 1 twist) for beam element corresponding to node 1
            matrix[12*element.node_index[0], 12*element.node_index[0]] = element.stretch_mod
            matrix[12*element.node_index[0]+1, 12*element.node_index[0]+1] = element.stretch_mod
            matrix[12*element.node_index[0]+2, 12*element.node_index[0]+2] = element.stretch_mod
            # 6 DoF (3 translational 2 bend 1 twist) for beam element corresponding to node 2
            matrix[12*element.node_index[1], 12*element.node_index[1]] = element.stretch_mod
            matrix[12*element.node_index[1]+1, 12*element.node_index[1]+1] = element.stretch_mod
            matrix[12*element.node_index[1]+2, 12*element.node_index[1]+2] = element.stretch_mod
            # Cross terms for beam element
            matrix[12*element.node_index[0], 12*element.node_index[1]] = element.stretch_mod
            matrix[12*element.node_index[0]+1, 12*element.node_index[1]+1] = element.stretch_mod
            matrix[12*element.node_index[0], 12*element.node_index[1]] = element.stretch_mod
            
    #if element is a torsional spring then there are 3 DoFs for rotations in x, y, z
        elif isinstance(element, Spring):
            # Node 1
            matrix[12*element.node_index[0]+3, 12*element.node_index[0]+3] = element.rot_x
            matrix[12*element.node_index[0]+4, 12*element.node_index[0]+4] = element.rot_y
            matrix[12*element.node_index[0]+5, 12*element.node_index[0]+5] = element.rot_z
            # Node 2
            matrix[12*element.node_index[1]+3, 12*element.node_index[1]+3] = element.rot_x
            matrix[12*element.node_index[1]+4, 12*element.node_index[1]+4] = element.rot_y
            matrix[12*element.node_index[1]+5, 12*element.node_index[1]+5] = element.rot_z
    #if element is a nick there will be reduced stretch, bend and twist moduli
    # single or double crossovers, open nicks and bulges modeled as torsional springs
    
    print(dx1dx1_e1, 'dx1dx1_e1 test')
    #print(pd.DataFrame(matrix))
    return matrix
 
test_system = [test_node_1, test_node_2, test_beam_1] 
 
 
 
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
 
generate_matrix(test_system)  
def main():
    matrix = generate_matrix(system)
    #visualise_spring_system(system)
    visualise(matrix, show=True)
    return
 
if __name__ == '__main__':
    pass
    #main()

