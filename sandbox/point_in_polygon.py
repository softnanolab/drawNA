import numpy as np
import matplotlib.pyplot as plt
from shapely import geometry

from drawNA.polygons import BoundaryPolygon

COLOURS = [
    (1, 0, 0, 1),
    (0, 1, 0, 1)
]

def square(L : int) -> np.ndarray:
    return np.array([
        [0., 0.],
        [1., 0.],
        [1., 1.],
        [0., 1.],
    ]) * L

def pentagon(L : int) -> np.ndarray:
    return (np.array([
        [0.550, 0.450],
        [0.455, 0.519],
        [0.491, 0.631],
        [0.609, 0.631],
        [0.645, 0.519]
    ]) - 0.5) * L

class PointInPolygon2D:
    def __init__(self, polygon_array, lattice, gridsize=0.1):
        self.polygon = np.concatenate([polygon_array, np.array([polygon_array[0, :]])], axis=0)
        self.lattice = lattice
        self.obj = BoundaryPolygon(polygon_array)
        _ = self.inside

    @property
    def box(self):
        return np.array([
            self.polygon[:, 0].min(),
            self.polygon[:, 0].max(),
            self.polygon[:, 1].min(),
            self.polygon[:, 1].max()
        ])

    def in_box(self, point):
        box = self.box
        return (
            point[0] > box[0] and \
            point[0] < box[1] and \
            point[1] > box[2] and \
            point[1] < box[3]
        )

    @property
    def inside(self):
        inside = []
        new_lattice =  []
        for i, point in enumerate(self.lattice):
            if not self.in_box(point):
                inside.append(False)
                continue

            polygon = geometry.Polygon(self.polygon)
            inside.append(polygon.contains(geometry.Point(point)))
            new_lattice.append(point)
        self.lattice = new_lattice
        return inside
        

    def plot(self):
        fig, ax = plt.subplots()
        box = np.array([
            [self.box[0], self.box[2]],
            [self.box[1], self.box[2]],
            [self.box[1], self.box[3]],
            [self.box[0], self.box[3]],
            [self.box[0], self.box[2]]
        ])
        ax.plot(self.polygon[:, 0], self.polygon[:, 1], 'k-')
        ax.plot(box[:, 0], box[:, 1], 'k:')
        for i, point in enumerate(self.lattice):
            ax.plot(
                point[0], 
                point[1],
                f'o', 
                markerfacecolor=COLOURS[int(self.inside[i])],
                markeredgecolor=COLOURS[int(self.inside[i])],
                )
        plt.show()

def generate_lattice(image_shape, lattice_vectors) :
    center_pix = np.array(image_shape) // 2
    # Get the lower limit on the cell size.
    dx_cell = max(abs(lattice_vectors[0][0]), abs(lattice_vectors[1][0]))
    dy_cell = max(abs(lattice_vectors[0][1]), abs(lattice_vectors[1][1]))
    # Get an over estimate of how many cells across and up.
    nx = image_shape[0]//dx_cell
    ny = image_shape[1]//dy_cell
    # Generate a square lattice, with too many points.
    # Here I generate a factor of 4 more points than I need, which ensures 
    # coverage for highly sheared lattices.  If your lattice is not highly
    # sheared, than you can generate fewer points.
    x_sq = np.arange(-nx, nx, dtype=float)
    y_sq = np.arange(-ny, nx, dtype=float)
    x_sq.shape = x_sq.shape + (1,)
    y_sq.shape = (1,) + y_sq.shape
    # Now shear the whole thing using the lattice vectors
    x_lattice = lattice_vectors[0][0]*x_sq + lattice_vectors[1][0]*y_sq
    y_lattice = lattice_vectors[0][1]*x_sq + lattice_vectors[1][1]*y_sq
    # Trim to fit in box.
    mask = ((x_lattice < image_shape[0]/2.0)
             & (x_lattice > -image_shape[0]/2.0))
    mask = mask & ((y_lattice < image_shape[1]/2.0)
                    & (y_lattice > -image_shape[1]/2.0))
    x_lattice = x_lattice[mask]
    y_lattice = y_lattice[mask]
    # Translate to the centre pix.
    x_lattice += center_pix[0]
    y_lattice += center_pix[1]
    # Make output compatible with original version.
    out = np.empty((len(x_lattice), 2), dtype=float)
    out[:, 0] = y_lattice
    out[:, 1] = x_lattice
    return out

def main(**kwargs):
    lattice = generate_lattice(
        (4, 4), 
        np.array([
            [0.1, 0.],
            [0., 0.1]
        ])) - 2.
    obj = PointInPolygon2D(
        pentagon(5), 
        lattice)
    obj.plot()
    return 

if __name__ == '__main__':
    main()