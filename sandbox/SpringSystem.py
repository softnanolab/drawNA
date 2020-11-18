import numpy as np
'''
The Node class creates a node object that has a position accessed by self.get_position. The set_position method allows the user
to input a numpy array with elements 0,1,2 representing the e0, e1 and e2 vectors. e0 is the position vector of the node, e1
points towards the major groove, e2 points towards the base in a chosen strand and e3 is automatically calculated by taking the 
vector product of e1 and e2
'''

class Node:

    def __init__(self):
        self.position = np.array([[0,0,0], [0,0,0], [0,0,0], [0,0,0]])
    def set_position(self, e0, e1, e2):
        '''
        Parameters: 
        Input: e0, e1, e2 (ndarray)  
        Output: [e0, e1, e2, e3] (ndarray) dim=4x3 
        '''
        self.position[0] = np.asarray(e0)
        self.position[1] = np.asarray(e1)
        self.position[2] = np.asarray(e2)
        self.position[3] = np.cross(self.position[1], self.position[2])
    def get_position(self):
        return self.positon

test_node = Node()
test_node.set_position([1,1,1], [0, -1, -1], [0, -2, -1])
print(type(test_node.position)) # returns ndarray
print(test_node.position)