#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "grid1D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##

"""
1D Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools.decorators import getsetDeprecated
from mesh1D import Mesh1D
from fipy.tools import parallel

from fipy.meshes.builders import Grid1DBuilder

class Grid1D(Mesh1D):
    """
    Creates a 1D grid mesh.
    
        >>> mesh = Grid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]
         
        >>> mesh = Grid1D(dx = (1, 2, 3))
        >>> print mesh.cellCenters
        [[ 0.5  2.   4.5]]
         
        >>> mesh = Grid1D(nx = 2, dx = (1, 2, 3))
        Traceback (most recent call last):
        ...
        IndexError: nx != len(dx)

    """
    def __init__(self, dx=1., nx=None, overlap=2, 
                 communicator=parallel,
                 BuilderClass=Grid1DBuilder):

        builder = BuilderClass()

        self.args = {
            'dx': dx, 
            'nx': nx, 
            'overlap': overlap
        }

        self.dim = 1
        self.dx = PhysicalField(value=dx)
        scale   = PhysicalField(value=1, unit=self.dx.unit)
        self.dx /= scale

        nx = self._calcNumPts(d=self.dx, n=nx)

        self.globalNumberOfCells = nx
        self.globalNumberOfFaces = nx + 1

        builder.buildParallelInfo([nx], overlap, communicator)

        ([self.nx],
         self.overlap,
         self.offset,
         self.occupiedNodes) = builder.getParallelInfo()

        if numerix.getShape(self.dx) is not ():
            Xoffset = numerix.sum(self.dx[0:self.offset])
            self.dx = self.dx[self.offset:self.offset + self.nx]
        else:
            Xoffset = self.dx * self.offset
            
        vertices              = self._createVertices() + ((Xoffset,),)
        self.numberOfVertices = len(vertices[0])
        faces                 = self._createFaces()
        self.numberOfFaces    = len(faces[0])
        cells                 = self._createCells()

        Mesh1D.__init__(self, vertices, faces, cells, communicator=communicator)
        
        self.scale = scale

    def __repr__(self):
        if self.args["nx"] is None:
            return "%s(dx=%s)" % (self.__class__.__name__, str(self.args["dx"]))
        else:
            return "%s(dx=%s, nx=%d)" % (self.__class__.__name__, str(self.args["dx"]), self.args["nx"])

    def _createVertices(self):
        x = self._calcVertexCoordinates(self.dx, self.nx)
        return x[numerix.newaxis,...]
    
    def _createFaces(self):
        if self.numberOfVertices == 1:
            return numerix.arange(0)[numerix.newaxis, ...]
        else:
            return numerix.arange(self.numberOfVertices)[numerix.newaxis, ...]

    def _createCells(self):
        """
        cells = (f1, f2) going left to right.
        f1 etc. refer to the faces
        """
        f1 = numerix.arange(self.nx)
        f2 = f1 + 1
        return numerix.array((f1, f2))

    @getsetDeprecated
    def getDim(self):
        return self.dim

    @getsetDeprecated
    def getScale(self):
        return self.scale['length']

    @getsetDeprecated
    def getPhysicalShape(self):
        return self.physicalShape

    @property
    def physicalShape(self):
        """Return physical dimensions of Grid1D.
        """
        from fipy.tools.dimensions.physicalField import PhysicalField
        return PhysicalField(value = (self.nx * self.dx * self.scale,))

    @getsetDeprecated
    def _getMeshSpacing(self):
        return self._meshSpacing

    @property
    def _meshSpacing(self):
        return numerix.array((self.dx,))[...,numerix.newaxis]
    
    @getsetDeprecated
    def getShape(self):
        return self.shape

    @property
    def shape(self):
        return (self.nx,)
        
    @getsetDeprecated
    def _getGlobalNonOverlappingCellIDs(self):
        return self._globalNonOverlappingCellIDs

    @property
    def _globalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1] for mesh A

            A        B
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """

        return numerix.arange(self.offset + self.overlap['left'], 
                              self.offset + self.nx - self.overlap['right'])

    @getsetDeprecated
    def _getGlobalOverlappingCellIDs(self):
        return self._globalOverlappingCellIDs

    @property
    def _globalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2] for mesh A

            A        B
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.offset, self.offset + self.nx)

    @getsetDeprecated
    def _getLocalNonOverlappingCellIDs(self):
        return self._localNonOverlappingCellIDs

    @property
    def _localNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1] for mesh A

            A        B
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.overlap['left'], 
                              self.nx - self.overlap['right'])

    @getsetDeprecated
    def _getLocalOverlappingCellIDs(self):
        return self._localOverlappingCellIDs

    @property
    def _localOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2] for mesh A

            A        B
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.nx)

    @getsetDeprecated
    def _getGlobalNonOverlappingFaceIDs(self):
        return self._globalNonOverlappingFaceIDs

    @property
    def _globalNonOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 2] for mesh A

            A    ||   B
        ------------------
        0   1    2   3   4
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.offset + self.overlap['left'], 
                              self.offset + self.numberOfFaces - self.overlap['right'])

    @getsetDeprecated
    def _getGlobalOverlappingFaceIDs(self):
        return self._globalOverlappingFaceIDs

    @property
    def _globalOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A    ||   B
        ------------------
        0   1    2   3   4
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.offset, self.offset + self.numberOfFaces)

    @getsetDeprecated
    def _getLocalNonOverlappingFaceIDs(self):
        return self._localNonOverlappingFaceIDs

    @property
    def _localNonOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 2] for mesh A

            A    ||   B
        ------------------
        0   1   2/1  2   3
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.overlap['left'], 
                              self.numberOfFaces - self.overlap['right'])

    @getsetDeprecated
    def _getLocalOverlappingFaceIDs(self):
        return self._localOverlappingFaceIDs

    @property
    def _localOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A   ||   B
        ------------------
        0   1   2   3    |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.numberOfFaces)

    
## pickling

    def __getstate__(self):
        return self.args
        
    def __setstate__(self, dict):
        self.__init__(**dict)

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. Fixed a bug where the following throws
        an error on solve() when nx is a float.

            >>> # from fipy import *
            >>> from fipy import CellVariable, DiffusionTerm
            >>> mesh = Grid1D(nx=3., dx=(1., 2., 3.))
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
