#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "uniformGrid1D.py"
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

from fipy.tools.numerix import MA
from fipy.meshes.topologies import _UniformMeshTopology1D
from fipy.meshes.geometries import _UniformGridGeometry1D

from grid1D import Grid1D
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix
from fipy.tools import parallel
from fipy.tools.decorators import getsetDeprecated

from fipy.variables.cellVariable import CellVariable
from fipy.variables.faceVariable import FaceVariable

class UniformGrid1D(Grid1D):
    """
    Creates a 1D grid mesh.
    
        >>> mesh = UniformGrid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]
         
    """
    def __init__(self, dx=1., nx=1, origin=(0,), overlap=2,
                       communicator=parallel,
                       GeomClass=_UniformGridGeometry1D):
        origin = numerix.array(origin)
        
        self.args = {
            'dx': dx, 
            'nx': nx, 
            'origin': origin, 
            'overlap': overlap
        }
        
        self.dim = 1
        
        self.dx = PhysicalField(value=dx)
        scale = PhysicalField(value=1, unit=self.dx.unit)
        self.dx /= scale
        
        nx = int(nx)

        self.globalNumberOfCells = nx
        self.globalNumberOfFaces = nx + 1
        
        (self.nx,
         self.overlap,
         self.offset) = self._calcParallelGridInfo(nx, overlap, communicator)
        
        self.origin = PhysicalField(value=origin)
        self.origin /= scale
        self.origin += self.offset * self.dx
        
        self.numberOfVertices = self.nx + 1
        if self.nx == 0:
            self.numberOfFaces = 0
        else:
            self.numberOfFaces = self.nx + 1
        self.numberOfCells = self.nx

        self._scale = {
            'length': 1.,
            'area': 1.,
            'volume': 1.
        }

        self._geometry = GeomClass(mesh=self,
                                   origin=self.origin,
                                   dx=self.dx,
                                   numberOfFaces=self.numberOfFaces,
                                   numberOfCells=self.numberOfCells,
                                   scale=self._scale)

        self._topology = _UniformMeshTopology1D(facesLeft=self.facesLeft, 
                                                facesRight=self.facesRight, 
                                                numberOfCells=self.numberOfCells, 
                                                numberOfFaces=self.numberOfFaces, 
                                                mesh=self)
        
        self.communicator = communicator

    def setScale(self, scale):
        self._setScale(scale)

    def _setScale(self, scale):
        self._geometry.scale = scale

    scale = property(lambda s: s._geometry.scale, _setScale)
        
    @getsetDeprecated
    def _getFaceAreas(self):
        return self._faceAreas

    def _translate(self, vector):
        return UniformGrid1D(dx=self.dx, 
                             nx=self.args['nx'], 
                             origin=self.args['origin'] + numerix.array(vector),
                             overlap=self.args['overlap'])

    def __mul__(self, factor):
        return UniformGrid1D(dx=self.dx * factor,
                             nx=self.args['nx'],
                             origin=self.args['origin'] * factor,
                             overlap=self.args['overlap'])

    def _getConcatenableMesh(self):
        from mesh1D import Mesh1D
        return Mesh1D(vertexCoords = self.vertexCoords, 
                      faceVertexIDs = self._createFaces(), 
                      cellFaceIDs = self._createCells())
                      
##     get topology methods

##         from common/mesh
        
    @getsetDeprecated
    def _getCellFaceIDs(self):
        return self.cellFaceIDs

    @property
    def cellFaceIDs(self):
        return CellVariable(mesh=self, 
                            value=MA.array(self._createCells()), 
                            elementshape=(2,))

    @property
    def _maxFacesPerCell(self):
        return 2
        
##         from numMesh/mesh

    @getsetDeprecated
    def _getVertexCoords(self):
        return self.vertexCoords

    @property
    def vertexCoords(self):
        return self.faceCenters

    @getsetDeprecated
    def getFaceCellIDs(self):
        return self.faceCellIDs

    @property
    def faceCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = MA.array((c1 - 1, c1))
        if self.numberOfFaces > 0:
            ids[0,0] = ids[1,0]
            ids[1,0] = MA.masked
            ids[1,-1] = MA.masked
        return FaceVariable(mesh=self, value=ids, elementshape=(2,))

    @getsetDeprecated
    def _getCellVertexIDs(self):
        return self._cellVertexIDs

    @property
    def _cellVertexIDs(self):
        c1 = numerix.arange(self.numberOfCells)
        return CellVariable(mesh=self, value=(c1 + 1, c1), elementshape=(2,))

##     scaling
    
    def _setScaledGeometry(self):
        pass

    def _getNearestCellID(self, points):
        """
        Test cases

           >>> from fipy import *
           >>> m = Grid1D(nx=3)
           >>> print m._getNearestCellID(([0., .9, 3.],))
           [0 0 2]
           >>> print m._getNearestCellID(([1.1],))
           [1]
           >>> m0 = Grid1D(nx=2, dx=1.)
           >>> m1 = Grid1D(nx=4, dx=.5)
           >>> print m0._getNearestCellID(m1.cellCenters.globalValue)
           [0 0 1 1]
           
        """
        nx = self.globalNumberOfCells
        
        if nx == 0:
            return numerix.arange(0)
            
        x0, = self.cellCenters.globalValue[...,0]        
        xi, = points
        dx = self.dx
        
        i = numerix.array(numerix.rint(((xi - x0) / dx)), 'l')
        i[i < 0] = 0
        i[i > nx - 1] = nx - 1

        return i

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. The following was broken, now fixed.

            >>> from fipy import *
            >>> mesh = Grid1D(nx=3., dx=1.)
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
