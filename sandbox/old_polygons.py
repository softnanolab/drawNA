
# ---plot---#
def plot_the_vertices(vertices: np.ndarray):
    fig = plt.figure()
    nodes = vertices
    triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1])
    # hello = triangulation.get_masked_triangles()
    plt.triplot(triangulation, "-k")
    plt.axes().set_aspect("equal", "datalim")
    fig.show()

    fig2 = plt.figure()
    new_nodes_x = np.append(nodes[:, 0], nodes[0, 0])
    new_nodes_y = np.append(nodes[:, 1], nodes[0, 1])
    print(new_nodes_x, new_nodes_y)
    plt.plot(new_nodes_x, new_nodes_y)
    plt.axes().set_aspect("equal", "datalim")
    fig2.show()

if __name__ == '__main__':
    print("old_polygons.py contains some functions from drawNA.polygons.")


# ---run the script---#
def make_polygon(polygon_vertices: np.ndarray, SHAPE=None):

    plot_the_vertices(polygon_vertices)

    i = 1
    if not SHAPE == None:
        # Current Functionalities
        print(f"This shape has {SHAPE.edges.__len__()} edges \n")
        print(
            f"The coordinates of edge number {i+1} are: {str(SHAPE.edges[i].vertices[0])} and {str(SHAPE.edges[i].vertices[1])} \n"
        )
        print(f"It has a size of {round(SHAPE.edges[i].length,4)}\n")
        print(f"This edge is a type of {SHAPE.edges[i].kind} edge \n")
        print(f"It's vector is {SHAPE.edges[i].vector}\n")
        print(f"It's unit vector is {SHAPE.edges[i].unit_vector}")
    return


def define_edges(vertices_of_polygon, index_1, index_2, edge_type):
    return Edge(vertices_of_polygon[index_1, :], vertices_of_polygon[index_2, :],)