#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.meshes.topologies import _MeshTopology
from fipy.meshes.geometries import _MeshGeometry
from fipy.meshes.abstractMesh import AbstractMesh

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import serial

class MeshAdditionError(Exception):
    pass
    
class Mesh(AbstractMesh):
    """Generic mesh class using numerix to do the calculations

        Meshes contain cells, faces, and vertices.

        This is built for a non-mixed element mesh.
    """

    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs, communicator=serial):
        """faceVertexIds and cellFacesIds must be padded with minus ones."""
         
        self.vertexCoords = vertexCoords
        self.faceVertexIDs = MA.masked_values(faceVertexIDs, -1)
        self.cellFaceIDs = MA.masked_values(cellFaceIDs, -1)
        self.communicator = communicator

        self.dim = self.vertexCoords.shape[0]

        if not hasattr(self, "numberOfFaces"):
            self.numberOfFaces = self.faceVertexIDs.shape[-1]
        if not hasattr(self, "numberOfCells"):
            self.numberOfCells = self.cellFaceIDs.shape[-1]
        if not hasattr(self, "globalNumberOfCells"):
            self.globalNumberOfCells = self.numberOfCells
        if not hasattr(self, "globalNumberOfFaces"):
            self.globalNumberOfFaces = self.numberOfFaces

        self.faceCellIDs = self._calcFaceCellIDs() 

        self._setTopology()
        self._setGeometry(scaleLength = 1.)


    """Topology methods"""
    def _setTopology(self):
        self._topology = _MeshTopology(self.cellFaceIDs, 
                                       self.faceCellIDs, 
                                       self.numberOfCells,
                                       self._maxFacesPerCell,
                                       self) # `self` only for int/ext face calc

    def _setGeometry(self, scaleLength = 1.):
        self._geometry = _MeshGeometry(self.dim,
                                       self.faceVertexIDs,
                                       self.vertexCoords,
                                       self.faceCellIDs,
                                       self.cellFaceIDs,
                                       self.numberOfCells,
                                       self._maxFacesPerCell,
                                       self._cellToFaceOrientations,
                                       scaleLength)
    
    @property
    def _concatenatedClass(self):
        return Mesh

    def __mul__(self, factor):
        """
        Dilate a `Mesh` by `factor`.
        
            >>> from fipy.meshes import Grid2D
            >>> baseMesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
            >>> print baseMesh.cellCenters
            [[ 0.5  1.5  0.5  1.5]
             [ 0.5  0.5  1.5  1.5]]

        The `factor` can be a scalar
        
            >>> dilatedMesh = baseMesh * 3
            >>> print dilatedMesh.cellCenters
            [[ 1.5  4.5  1.5  4.5]
             [ 1.5  1.5  4.5  4.5]]

        or a vector
        
            >>> dilatedMesh = baseMesh * ((3,), (2,))
            >>> print dilatedMesh.cellCenters
            [[ 1.5  4.5  1.5  4.5]
             [ 1.   1.   3.   3. ]]

        
        but the vector must have the same dimensionality as the `Mesh`
        
            >>> dilatedMesh = baseMesh * ((3,), (2,), (1,))
            Traceback (most recent call last):
            ...
            ValueError: shape mismatch: objects cannot be broadcast to a single shape
        """ 
        newCoords = self.vertexCoords * factor
        newmesh = Mesh(vertexCoords=newCoords, 
                       faceVertexIDs=numerix.array(self.faceVertexIDs), 
                       cellFaceIDs=numerix.array(self.cellFaceIDs))
        return newmesh

    __rmul__ = __mul__
    
    @property
    def _concatenableMesh(self):
        return self

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    """calc Topology methods"""

    """
    TODO: Does this really belong in mesh? I don't think so.
    """
    def _calcFaceCellIDs(self):
        array = MA.array(MA.indices(self.cellFaceIDs.shape, 'l')[1], 
                         mask=MA.getmask(self.cellFaceIDs))
        faceCellIDs = MA.zeros((2, self.numberOfFaces), 'l')

        ## Nasty bug: MA.put(arr, ids, values) fills its ids and
        ## values arguments when masked!  This was not the behavior
        ## that was assumed when used below.  It was only working
        ## because the old fill value was 0 and the first element of
        ## the array needed to be 0 since the cell's face was
        ## 0. numerix.put() has been changed to deal with this
        ## properly.

##         MA.put(firstRow, cellFaceIDsFlat[::-1], array[::-1])
##         MA.put(secondRow, cellFaceIDsFlat, array)
        firstRow = faceCellIDs[0]
        secondRow = faceCellIDs[1]
        numerix.put(firstRow, self.cellFaceIDs[::-1,::-1], array[::-1,::-1])
        numerix.put(secondRow, self.cellFaceIDs, array)
        
        mask = ((False,) * self.numberOfFaces, (firstRow == secondRow))
        return MA.sort(MA.array(faceCellIDs, mask = mask),
                                   axis=0)

    """get Topology methods"""
    
    @property
    def _maxFacesPerCell(self):
        return self.cellFaceIDs.shape[0]

    def _isOrthogonal(self):
        return False
       
    @property
    def _cellVertexIDs(self):
        ## Get all the vertices from all the faces for each cell
        cellFaceVertices = numerix.take(self.faceVertexIDs, self.cellFaceIDs, axis=1)

        ## get a sorted list of vertices for each cell 
        cellVertexIDs = numerix.reshape(cellFaceVertices, (-1, self.numberOfCells))
        cellVertexIDs = MA.sort(cellVertexIDs, axis=0, fill_value=-1)

        cellVertexIDs = MA.sort(MA.concatenate((cellVertexIDs[-1, numerix.newaxis], 
                                                MA.masked_where(cellVertexIDs[:-1] 
                                                                == cellVertexIDs[1:], 
                                                                cellVertexIDs[:-1]))), 
                                axis=0, fill_value=-1)
        
        ## resize the array to remove extra masked values
        if cellVertexIDs.shape[-1] == 0:
            length = 0
        else:
            length = min(numerix.sum(MA.getmaskarray(cellVertexIDs), axis=0))
        return cellVertexIDs[length:][::-1]

    """
    Below is an ordered version of _getCellVertexIDs()
    It works for the test case in this file (other than the ordering, obviously)
    I've left it in as it may be useful when we need ordered vertices for cells
  
    def _getOrderedCellVertexIDs(self):

        ## Get all the vertices from all the faces for each cell
        from fipy.tools.numerix import take
        cellFaceVertices = take(self.faceVertexIDs, self.cellFaceIDs)

        ## get a sorted list of vertices for each cell
        NCells = self.numberOfCells
        cellVertexIDs = MA.reshape(cellFaceVertices.flat, (NCells, -1))
        newmask = MA.getmaskarray(cellVertexIDs).copy()

        for i in range(len(cellVertexIDs[0]) - 1):
            for j in range(len(cellVertexIDs[0]))[i + 1:]:

                newmask[:,j] = MA.where(newmask[:,j],
                                        newmask[:,j],
                                        MA.where(MA.filled(cellVertexIDs)[:,i] == MA.filled(cellVertexIDs)[:,j],
                                                 1,
                                                 newmask[:,j]))


        cellVertexIDs = MA.masked_array(cellVertexIDs, newmask)

        for i in range(len(cellVertexIDs[0]) - 1):
            j = i + 1
            while j < len(cellVertexIDs[0]):
                tmp = cellVertexIDs[:]
                tmp[:, i] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
                                     MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
                                              cellVertexIDs[:, i],
                                              cellVertexIDs[:, j]),
                                     cellVertexIDs[:, i])
                                                      
                tmp[:, j] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
                                     MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
                                              cellVertexIDs[:,j],
                                              cellVertexIDs[:,i]),
                                     cellVertexIDs[:, j])

                cellVertexIDs = tmp[:]

                j += 1


        length = len(cellVertexIDs[0]) - min(numerix.sum(MA.getmaskarray(cellVertexIDs), axis = 1))
        return cellVertexIDs[:, :length]
    """


    """scaling"""
    
    def _getNearestCellID(self, points):
        """
        Test cases

           >>> from fipy import *
           >>> m0 = Grid2D(dx=(.1, 1., 10.), dy=(.1, 1., 10.))
           >>> m1 = Grid2D(nx=2, ny=2, dx=5., dy=5.)
           >>> print m0._getNearestCellID(m1.cellCenters.globalValue)
           [4 5 7 8]
           
        """
        return numerix.nearest(data=self.cellCenters.globalValue, points=points)
        

    """pickling"""

    def __getstate__(self):
        dict = {
            'vertexCoords' : self.vertexCoords *  self.scale['length'],            
            'faceVertexIDs' : self.faceVertexIDs,
            'cellFaceIDs' : self.cellFaceIDs }
        return dict

    def __setstate__(self, dict):
        Mesh.__init__(self, **dict)
##        self.__init__(dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])
     
    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
            >>> dx = 2.
            >>> dy = 1.23456
            >>> dz = 1.e-1
            
            >>> vertices = numerix.array(((0., 1., 1., 0., 0., 1., 1., 0., 2., 2.),
            ...                           (0., 0., 1., 1., 0., 0., 1., 1., 0., 0.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 0., 1.)))
            >>> vertices *= numerix.array([[dx], [dy], [dz]])
            >>> faces = MA.masked_values(((0, 7, 3, 5, 1, 3, 1, 9, 8, 8),
            ...                           (1, 6, 7, 6, 0, 2, 8, 5, 1, 9),
            ...                           (2, 5, 4, 2, 4, 6, 2, 6, 5, 6),
            ...                           (3, 4, 0, 1, 5, 7, -1, -1, 9, 2)), -1)
            >>> cells = MA.masked_values(((0, 3),
            ...                           (1, 6), 
            ...                           (2, 7),
            ...                           (3, 8),
            ...                           (4, 9),
            ...                           (5, -1)), -1)

            >>> mesh = Mesh(vertexCoords=vertices, faceVertexIDs=faces, cellFaceIDs=cells)

            >>> externalFaces = numerix.array((0, 1, 2, 4, 5, 6, 7, 8, 9))
            >>> print numerix.allequal(externalFaces, 
            ...                        numerix.nonzero(mesh.exteriorFaces))
            1

            >>> internalFaces = numerix.array((3,))
            >>> print numerix.allequal(internalFaces, 
            ...                        numerix.nonzero(mesh.interiorFaces))
            1

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0,  0,  0, 0,  0,  0,  1,  1,  1,  1),
            ...                                 (-1, -1, -1, 1, -1, -1, -1, -1, -1, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.faceCellIDs)
            1
            
            >>> dxdy = dx * dy
            >>> dxdz = dx * dz
            >>> dydz = dy * dz
            >>> faceAreas = numerix.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
            ...                            dxdy/2., dxdy/2., dxdz, numerix.sqrt(dx**2 + dy**2) * dz))
            >>> numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, MA.filled(faces, 0), axis=1)
            >>> faceCenters = faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]
            >>> numVex = numerix.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
            >>> faceCenters /= numVex
            >>> numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array((( 0., 0., -1., 1.,  0., 0.,  0., 0.,  0., dy / numerix.sqrt(dy**2 + dx**2)),
            ...                              ( 0., 0.,  0., 0., -1., 1.,  0., 0., -1., dx / numerix.sqrt(dy**2 + dx**2)),
            ...                              (-1., 1.,  0., 0.,  0., 0., -1., 1.,  0., 0.)))
            >>> numerix.allclose(faceNormals, mesh._faceNormals, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, -1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, -2)), -2)
            >>> numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations)
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz / 2.))
            >>> numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2., dx+dx/3.),
            ...                              (dy/2.,    dy/3.),
            ...                              (dz/2.,    dz/2.)))
            >>> print numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> d1 = numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2)
            >>> d2 = numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2)
            >>> d3 = numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)
            >>> d4 = numerix.sqrt((5 * dx / 6.)**2 + (dy / 6.)**2)
            >>> faceToCellDistances = MA.masked_values(((dz / 2., dz / 2., dx / 2., dx / 2., dy / 2., dy / 2., dz / 2., dz / 2., d2, d3),
            ...                                         (     -1,      -1,      -1,      d1,      -1,      -1,      -1,      -1, -1, -1)), -1)
            >>> print numerix.allclose(faceToCellDistances, mesh._faceToCellDistances, atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> cellDistances = numerix.array((dz / 2., dz / 2., dx / 2.,
            ...                                d4,
            ...                                dy / 2., dy / 2., dz / 2., dz / 2.,
            ...                                d2,
            ...                                d3))
            >>> numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._faceToCellDistanceRatio, atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> numerix.allclose(areaProjections, mesh._areaProjections, atol = 1e-10, rtol = 1e-10)
            1

            >>> v1 = numerix.take(vertices, numerix.array(faces[0]), axis=1)
            >>> tmp = faceCenters - v1
            >>> tangents1 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> numerix.allclose(tangents1, mesh._faceTangents1, atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp = numerix.cross(tangents1, faceNormals, axis=0)
            >>> tangents2 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> numerix.allclose(tangents2, mesh._faceTangents2, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1,  0),
            ...                                   (-1, -1),
            ...                                   (-1, -1), 
            ...                                   ( 1, -1), 
            ...                                   (-1, -1),
            ...                                   (-1, -1)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._cellToCellIDs)
            1

            >>> cellToCellDistances = MA.masked_values(((dz / 2., d4),
            ...                                         (dz / 2., -1),
            ...                                         (dx / 2., -1),
            ...                                         (     d4, -1),
            ...                                         (dy / 2., -1),
            ...                                         (dy / 2., -1)), -1)
            >>> numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._interiorCellIDs)
            1

            >>> exteriorCellIDs = numerix.array((0, 1))
            >>> numerix.allequal(exteriorCellIDs, mesh._exteriorCellIDs)
            1

            >>> cellNormals = MA.masked_values((((0, -1),
            ...                                  (0, 0),
            ...                                  (-1, 0),
            ...                                  (1, 0),
            ...                                  (0, dy / numerix.sqrt(dx**2 + dy**2)),
            ...                                  (0, -1000)),
            ...                                 ((0, 0),
            ...                                  (0, 0),
            ...                                  (0, 0),
            ...                                  (0, -1),
            ...                                  (-1, dx / numerix.sqrt(dx**2 + dy**2)),
            ...                                  (1, -1000)),
            ...                                 ((-1, 0),
            ...                                  (1, -1),
            ...                                  (0, 1),
            ...                                  (0, 0),
            ...                                  (0, 0),
            ...                                  (0, -1000))), -1000)
            >>> numerix.allclose(cellNormals, mesh._cellNormals, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellAreaProjections = MA.masked_values((((0, -dy * dz),
            ...                                          (0, 0),
            ...                                          (-dy * dz, 0),
            ...                                          ( dy * dz, 0),
            ...                                          (0, dy * dz),
            ...                                          (0, -1000)),
            ...                                         ((0, 0),
            ...                                          (0, 0),
            ...                                          (0, 0),
            ...                                          (0, -dx * dz),
            ...                                          (-dx * dz, dx * dz),
            ...                                          ( dx * dz, -1000)),
            ...                                         ((-dx * dy, 0),
            ...                                          ( dx * dy, -dx * dy / 2.),
            ...                                          (0, dx * dy / 2.),
            ...                                          (0, 0),
            ...                                          (0, 0),
            ...                                          (0, -1000))), -1000)
            >>> numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_values(((7, 6, 5, 4, 3, 2, 1, 0), (9, 8, 6, 5, 2, 1, -1000, -1000)), -1000)
            >>> cellVertexIDs = MA.masked_values(((7, 9),
            ...                                   (6, 8),
            ...                                   (5, 6),
            ...                                   (4, 5),
            ...                                   (3, 2),
            ...                                   (2, 1),
            ...                                   (1, -1000),
            ...                                   (0, -1000)), -1000)
            >>> numerix.allclose(cellVertexIDs, mesh._cellVertexIDs)
            1


            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allequal(mesh.cellCenters, unpickledMesh.cellCenters)
            True

            >>> dx = 1.
            >>> dy = 1.
            >>> nx = 10
            >>> ny = 2
            >>> from fipy.meshes import Grid2D
            >>> gridMesh = Grid2D(dx, dy, nx, ny)
            >>> from fipy.meshes.tri2D import Tri2D
            >>> triMesh = Tri2D(dx, dy, nx, 1) + [[dx*nx], [0]]
            >>> bigMesh = gridMesh + triMesh
            >>> x, y = bigMesh.cellCenters
            >>> from fipy.variables.cellVariable import CellVariable
            >>> volumes = CellVariable(mesh=bigMesh, value=1.)
            >>> volumes[x > dx * nx] = 0.25
            >>> print numerix.allclose(bigMesh.cellVolumes, volumes)
            True
            
            Following test was added due to a bug in adding UniformGrids.

            >>> from fipy.meshes.uniformGrid1D import UniformGrid1D
            >>> a = UniformGrid1D(nx=10) + (10,)
            >>> print a.cellCenters
            [[ 10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5]]
            >>> b = 10 + UniformGrid1D(nx=10)
            >>> print b.cellCenters
            [[ 10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5]]
            
            >>> from fipy.tools import parallel
            >>> if parallel.Nproc == 1:
            ...     c =  UniformGrid1D(nx=10) + (UniformGrid1D(nx=10) + 10)
            >>> print (parallel.Nproc > 1 
            ...        or numerix.allclose(c.cellCenters[0],
            ...                            [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
            ...                            12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]))
            True

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
