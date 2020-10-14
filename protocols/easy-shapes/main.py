import numpy as np
from drawNA.polygons import BoundaryPolygon, make_polygon

###---stored arrays of some shapes---###
square = np.array([[0, 0, 0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
trapezium = np.array([[0, 0, 0], [1.0, 10.0, 0.0], [11.0, 10.0, 0.0], [15.0, 0.0, 0.0]])

no = 0.866025
hexagon = np.array(
    [[0, 1, 0], [no, 0.5, 0], [no, -0.5, 0], [0, -1, 0], [-no, -0.5, 0], [-no, 0.5, 0]]
)

star = np.array(
    [
        [0, 300, 0],
        [100, 110, 0],
        [300, 70, 0],
        [160, -90, 0],
        [190, -300, 0],
        [0, -210, 0],
        [-190, -300, 0],
        [-160, -90, 0],
        [-300, 70, 0],
        [-100, 110, 0],
    ]
)


###---run the program---###
SHAPE = BoundaryPolygon(star)
make_polygon(star, SHAPE)
