import numpy as np

class Node:
    def __init__(self):
        self.position = np.array([[0,0,0], [0,0,0], [0,0,0], [0,0,0]])
    def set_position(self, *args):
        self.position[0] = np.asarray(args[0])
        self.position[1] = np.asarray(args[1])
        self.position[2] = np.asarray(args[2])
        self.position[3] = np.cross(self.position[1], self.position[2])
    def get_position(self):
        return self.positon

test_node = Node()
test_node.set_position([1,1,1], [0, -1, -1], [0, -2, -1])
print(type(test_node.position)) 
print(test_node.position)