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

from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix
from fipy.tools import parallel

from fipy.meshes.builders import UniformGrid1DBuilder
from fipy.meshes.builders import Grid1DBuilder
from fipy.meshes.abstractGrid import AbstractGrid1DFactory
from fipy.meshes.abstractMesh import AbstractMesh

class UniformGrid1D(AbstractGrid1DFactory(AbstractMesh)):
    """
    Creates a 1D grid mesh.
    
        >>> mesh = UniformGrid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]
         
    """
    def __init__(self, dx=1., nx=1, origin=(0,), overlap=2,
                       communicator=parallel,
                       GeomClass=_UniformGridGeometry1D):
        builder = UniformGrid1DBuilder()

        origin = numerix.array(origin)
        
        self.args = {
            'dx': dx, 
            'nx': nx, 
            'origin': origin, 
            'overlap': overlap
        }
        
        self._scale = {
            'length': 1.,
            'area': 1.,
            'volume': 1.
        }
         
        builder.buildGridData([dx], [nx], overlap, communicator, origin)

        ([self.dx],
         [self.nx],
         self.dim,
         scale,
         self.globalNumberOfCells,
         self.globalNumberOfFaces,
         self.overlap,
         self.offset,
         self.numberOfVertices,
         self.numberOfFaces,
         self.numberOfCells,
         self.shape,
         self.physicalShape,
         self._meshSpacing,
         self.occupiedNodes,
         self.origin) = builder.gridData

        self._setGeometry(GeomClass=GeomClass)
        self._setTopology()

        self.communicator = communicator

    def _setTopology(self):
        self._topology = _UniformMeshTopology1D(self.facesLeft, 
                                                self.facesRight, 
                                                self.numberOfCells, 
                                                self.numberOfFaces, 
                                                self)

    def _setGeometry(self, GeomClass=_UniformGridGeometry1D):
        self._geometry = GeomClass(self.origin,
                                   self.dx,
                                   self.numberOfFaces,
                                   self.numberOfCells,
                                   scale=self._scale)
                                                          
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

    @property
    def _concatenableMesh(self):
        from mesh1D import Mesh1D
        return Mesh1D(vertexCoords = self.vertexCoords, 
                      faceVertexIDs = Grid1DBuilder.createFaces(self.numberOfVertices), 
                      cellFaceIDs = Grid1DBuilder.createCells(self.nx))
                      
##     get topology methods

##         from common/mesh

    @property
    def _cellFaceIDs(self):
        return MA.array(Grid1DBuilder.createCells(self.nx))

    @property
    def _maxFacesPerCell(self):
        return 2
        
##         from numMesh/mesh

    @property
    def vertexCoords(self):
        return self.faceCenters

    @property
    def faceCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = MA.array((c1 - 1, c1))
        if self.numberOfFaces > 0:
            ids[0,0] = ids[1,0]
            ids[1,0] = MA.masked
            ids[1,-1] = MA.masked
        return ids

    @property
    def _cellVertexIDs(self):
        c1 = numerix.arange(self.numberOfCells)
        return numerix.array((c1 + 1, c1))

##     scaling
    
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
