from __future__ import unicode_literals
from builtins import object
from builtins import range
__all__ = []

class _Vertex(object):
    def __init__(self, ID, x, y):
        self.ID = ID
        self.up = None
        self.down = None
        self.closeVertices = None
        self.inLine = False
        self.x = x
        self.y = y

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def distance(self, vertex):
        from fipy.tools import numerix
        return numerix.sqrt((vertex.getX() - self.getX())**2 + (vertex.getY() - self.getY())**2)

    def getID(self):
        return self.ID

    def setCloseVertices(self, vertices):
        self.closeVertices = vertices

    def getCloseVertices(self):
        return self.closeVertices

    def getUp(self):
        return self.up

    def getDown(self):
        return self.down

    def setInLineTrue(self):
        self.inLine = True

    def getInLine(self):
        return self.inLine

    def setUp(self, vertex):
        self.setInLineTrue()
        self.up = vertex

    def setDown(self, vertex):
        self.setInLineTrue()
        self.down = vertex

class _Line(object):
    def __init__(self, seedVertex):
        if seedVertex.getUp() is not None or \
           seedVertex.getDown() is not None or \
           seedVertex.getInLine():
            raise ValueError('bad seedVertex')

        setVs = []

        vertex = seedVertex
        while vertex is not None and vertex.getUp() is None:
            for v in vertex.getCloseVertices():
                if v not in setVs:
                    if vertex in v.getCloseVertices():
                        if v.getDown() is None:
                            if vertex.getUp() is None:
                                if vertex.getDown() is not v:
                                    vertex.setUp(v)
                                    v.setDown(vertex)

            setVs.append(vertex)
            vertex = vertex.getUp()



        self.startVertex = seedVertex
        vertex = seedVertex
        while vertex is not None and vertex.getDown() is None:
            for v in vertex.getCloseVertices():
                if v not in setVs:
                    if vertex in v.getCloseVertices():
                        if v.getUp() is None:
                            if vertex.getDown() is None:
                                if vertex.getUp() is not v:
                                    vertex.setDown(v)
                                    v.setUp(vertex)

            self.startVertex = vertex
            setVs.append(vertex)
            vertex = vertex.getDown()

        self.startVertex.setInLineTrue()

    def getVertexListIDs(self):
        vertexListIDs = []
        vertexListIDs.append(self.startVertex.getID())
        v = self.startVertex.getUp()
        while v is not None and v is not self.startVertex:
            vertexListIDs.append(v.getID())
            v = v.getUp()

        if v is not None:
            vertexListIDs.append(v.getID())
        return vertexListIDs

def _getOrderedLines(IDs, coordinates, thresholdDistance = 0.0):
    """
    This function takes a set of IDs and corresponding coordinates and makes
    a set of closed curves.

    The following are a general set of test cases.

       >>> _getOrderedLines((0, 1, 2, 3), ((0, 0), (2, 0), (0, 1), (2, 1)))
       [[0, 2, 3, 1]]
       >>> _getOrderedLines((0, 1, 2, 3, 4), ((-10, -10), (0, 0), (2, 0), (0, 1), (2, 1)))
       [[0], [1, 3, 4, 2]]
       >>> _getOrderedLines((0, 1, 2, 3), ((0, 0), (0.9, 0), (0, 1), (1, 1)))
       [[0, 1, 3, 2]]
       >>> _getOrderedLines((0, 1, 2, 3, 4, 5), ((0, 0), (1, 0), (2, 0), (0, 1.1), (1, 1.1), (2, 1.1)))
       [[0, 1, 2, 5, 4, 3]]
       >>> _getOrderedLines((4, 5, 0, 1, 3, 2, 6), ((0, 0), (1, 0), (3, 0), (0, 1.1), (1, 1.1), (3, 1), (4, 1)))
       [[4, 5, 3, 1], [0, 2, 6]]
       >>> _getOrderedLines((4, 5, 3, 2, 1, 0, 9, 8, 7, 6), ((0, 0), (1, 0), (0, 1.1), (1, 1.1), (0, 3), (1, 3), (-1, 4), (2, 4), (-2, 4), (3, 4)))
       [[4, 5, 2, 3], [7, 9, 1, 0, 8, 6]]
       >>> from builtins import range
       >>> _getOrderedLines(list(range(7)), ((-7, 0), (-6, 0), (-5, 0), (0, 0), (5, 0), (6, 0), (7, 0)))
       [[0, 1, 2], [3], [4, 5, 6]]
       >>> from builtins import range
       >>> _getOrderedLines(list(range(7)), ((-7, 0), (-6, 0), (-5, 0), (0, 0), (5, 0), (6, 0), (7, 0)), thresholdDistance = 5.5)
       [[0, 1, 2, 3, 4, 5, 6]]

    Parameters
    ----------
    IDs : ndarray of int
    coordinates : ndarray
        An array of coordinates of the same length as IDs.
    """

    from fipy.tools import numerix
    coordinates = numerix.array(coordinates)
    closeIDs = numerix.zeros((len(IDs), len(IDs)), 'l')
    vertices = []
    for ID in IDs:
        distances = numerix.zeros(len(IDs), 'd')
        for i in range(len(coordinates[0,:])):
            distances += (coordinates[:, i] - coordinates[ID, i])**2
        vertices.append(_Vertex(ID, coordinates[ID, 0], coordinates[ID, 1]))
        closeIDs[ID,:] = numerix.argsort(distances)


    for ID in IDs:
        i = 1
        closeVertices = []
        while i < 3 or vertices[ID].distance(vertices[closeIDs[ID, i]]) < thresholdDistance:
            closeVertices.append(vertices[closeIDs[ID, i]])
            i += 1

        vertices[ID].setCloseVertices(closeVertices)

    listOfVertexLists = []

    for vertex in vertices:
        if not vertex.getInLine():
            listOfVertexLists.append(_Line(vertex).getVertexListIDs())

    return listOfVertexLists

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
