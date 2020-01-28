from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import random
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import serialComm

from fipy.meshes.mesh2D import Mesh2D
from fipy.meshes import Grid2D

__all__ = ["SkewedGrid2D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SkewedGrid2D(Mesh2D):
    """
    Creates a 2D grid mesh with horizontal faces numbered first and then
    vertical faces.  The points are skewed by a random amount (between `rand`
    and `-rand`) in the X and Y directions.

    .. note:: This `Mesh` only operates in serial
    """
    def __init__(self, dx = 1., dy = 1., nx = None, ny = 1, rand = 0, *args, **kwargs):
        self.args = {
            'dx': dx,
            'dy': dy,
            'nx': nx,
            'ny': ny,
            'rand': rand
        }

        self.nx = nx
        self.ny = ny

        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value = 1, unit = self.dx.unit)
        self.dx /= scale

        self.dy = PhysicalField(value = dy)
        if self.dy.unit.isDimensionless():
            self.dy = dy
        else:
            self.dy /= scale

        self.grid = Grid2D(nx=nx, ny=ny, dx=dx, dy=dy, communicator=serialComm)

        self.numberOfVertices = self.grid._numberOfVertices

        vertices = self.grid.vertexCoords

        changedVertices = numerix.zeros(vertices.shape, 'd')

        for i in range(len(vertices[0])):
            if((i % (nx+1)) != 0 and (i % (nx+1)) != nx and (i // nx + 1) != 0 and (i // nx + 1) != ny):
                changedVertices[0, i] = vertices[0, i] + (rand * ((random.random() * 2) - 1))
                changedVertices[1, i] = vertices[1, i] + (rand * ((random.random() * 2) - 1))
            else:
                changedVertices[0, i] = vertices[0, i]
                changedVertices[1, i] = vertices[1, i]


        faces = self.grid.faceVertexIDs

        cells = self.grid.cellFaceIDs

        Mesh2D.__init__(self, changedVertices, faces, cells, communicator=serialComm, *args, **kwargs)

        self.scale = scale

    @property
    def physicalShape(self):
        """Return physical dimensions of `Grid2D`.
        """
        return PhysicalField(value = (self.nx * self.dx * self.scale, self.ny * self.dy * self.scale))

    @property
    def _meshSpacing(self):
        return numerix.array((self.dx, self.dy))[..., numerix.newaxis]

    @property
    def shape(self):
        return (self.nx, self.ny)
