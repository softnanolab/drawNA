from typing import Tuple, Iterable

import numpy as np
import pandas as pd

import matplotlib.animation
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

class Node:
    def __init__(self, position, e1, e3):
        self.position = np.array(position)
        self._e1 = np.array(e1)
        self._e3 = np.array(e3)
        
    @property
    def _e2(self):
        return np.cross(self._e1, self._e3)
    
    def add_to_ax(self, ax: plt.Axes, **kwargs):
        point, = ax.plot(self.position[0], self.position[1], self.position[2], color='k', marker='o', markersize=10, **kwargs)
        arrows = []
        for vector, colour in zip([self._e1, self._e2, self._e3], ['r', 'g', 'b']):
            quiver = ax.quiver(*self.position, *vector, color=colour, length=0.5, **kwargs)
            arrows.append(quiver)
        return [point, arrows[0], arrows[1], arrows[2]]
    
class Connector:
    def __init__(self, nodes: Tuple[Node, Node]):
        self.nodes = nodes
        
    def add_to_ax(self, ax: plt.Axes):
        line, = ax.plot(
            [self.nodes[0].position[0], self.nodes[1].position[0]],
            [self.nodes[0].position[1], self.nodes[1].position[1]],
            [self.nodes[0].position[2], self.nodes[1].position[2]],
            'k:'
        )
        return line

class Animation:
    def __init__(self, nodes = None, connectors = None, limits=((-0.5, 2.5), (-0.5, 2.5), (-1.5, 1.0))):
        print('Creating animation')
        self.nodes = nodes
        self.connectors = connectors
        self.artists = []
        self.limits = limits
        self.setup()

    def set_limits(self):
        if not isinstance(self.limits, Iterable):
            limits = [self.limits, self.limits]
            self.ax.set_xlim(*limits)
            self.ax.set_ylim(*limits)
            self.ax.set_zlim(*limits)
        elif len(self.limits) == 1:
            limits = [self.limits] * 2
            self.ax.set_xlim(*limits)
            self.ax.set_ylim(*limits)
            self.ax.set_zlim(*limits)
        elif len(self.limits) == 2:
            self.ax.set_xlim(*self.limits)
            self.ax.set_ylim(*self.limits)
            self.ax.set_zlim(*self.limits)
        elif len(self.limits) == 3:
            self.ax.set_xlim(*self.limits[0])
            self.ax.set_ylim(*self.limits[1])
            self.ax.set_zlim(*self.limits[2])

    def draw(self):
        """initialize animation"""
        for i, node in enumerate(self.nodes):
            if i != 3:
                self.artists += node.add_to_ax(self.ax)
            else:
                self.artists += node.add_to_ax(self.ax, alpha=0.25)
        for connector in self.connectors:
            self.artists.append(connector.add_to_ax(self.ax))
        return self.artists
    
    def setup(self):
        self.artists = []
        self.figure = plt.figure()
        self.ax = self.figure.gca(projection='3d')
        self.ax.set_axis_off()
        self.ax.set_box_aspect(aspect = (1,1,1))
        self.set_limits()
        return
    
nodes = [
    Node(
        np.array([0., 1., 0.]),
        np.array([1., 0., 0.]),
        np.array([0., 0., 1.]),
    ),
    Node(
        np.array([1., 1., 0.]),
        np.array([1., 0., 0.]),
        np.array([0., 0., 1.]),
    ),
    Node(
        np.array([2., 1., 0.]),
        np.array([1., 0., 0.]),
        np.array([0., 0., 1.]),
    ),
    Node(
        np.array([1., 1., 0.]),
        np.array([1., 0., 0.]),
        np.array([0., 0., 1.]),
    ),
]

connectors = [
    Connector((nodes[:2])),
    Connector((nodes[1:])),
]

ani = Animation(nodes=nodes, connectors=connectors, limits=(0.5, 2.5))

def animate(t):    
    global ani
    ani.ax.cla()
    ani.ax.set_axis_off()
    ani.set_limits()
    ani.ax.set_box_aspect(aspect = (1,1,1))
    ani.nodes[1].position[1] += 0.25 * np.cos(t)
    return ani.draw()

if True:
    animation = matplotlib.animation.FuncAnimation(
        ani.figure, 
        animate, 
        frames=100, 
        interval=150, 
        init_func=ani.draw,
        repeat=False,
    )
    animation.save('animation.gif', writer='imagemagick', fps=10)

plt.show()
