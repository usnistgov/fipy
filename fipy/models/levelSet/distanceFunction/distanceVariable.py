#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "distanceVariable.py"
 #                                    created: 7/29/04 {10:39:23 AM} 
 #                                last update: 7/29/04 {11:03:01 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

import Numeric

from fipy.variables.cellVariable import CellVariable
import fipy.tools.array as array

class DistanceVariable(CellVariable):
    """

    The 'DistanceVariable` evaluates quantities associated with the
    distance function. It is mainly evaluated in the
    'DistanceFunctionEquation`.


    """
    
    def getCellInterfaceAreas(self):
        """
        Returns the length of the interface that crosses the cell

        A simple 1D test:

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 4, ny = 1)
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-1.5, -0.5, 0.5, 1.5))
           >>> Numeric.allclose(distanceVariable.getCellInterfaceAreas(), (0, 0., 1., 0))
           1

        A 2D test case:

           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (1.5, 0.5, 1.5,
           ...                                                          0.5,-0.5, 0.5,
           ...                                                          1.5, 0.5, 1.5))
           >>> Numeric.allclose(distanceVariable.getCellInterfaceAreas(), (0, 1, 0, 1, 0, 1, 0, 1, 0))
           1

        Another 2D test case:

           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-0.5, 0.5, 0.5, 1.5))
           >>> Numeric.allclose(distanceVariable.getCellInterfaceAreas(), (0, Numeric.sqrt(2) / 4,  Numeric.sqrt(2) / 4, 0))
           1

        Test to check that the circumfrence of a circle is in face 2*pi*r

           >>> mesh = Grid2D(dx = 0.05, dy = 0.05, nx = 20, ny = 20)
           >>> r = 0.25
           >>> rad = Numeric.sqrt((mesh.getCellCenters()[:,0] - .5)**2 + (mesh.getCellCenters()[:,1] - .5)**2) - r
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = rad)
           >>> print Numeric.sum(distanceVariable.getCellInterfaceAreas())
           1.57984690073
           
        """

        return Numeric.sum(abs(array.dot(self.getCellInterfaceNormals(), self.mesh.getCellAreaProjections(), axis = 2)), axis = 1)

    def getCellInterfaceNormals(self):
        """
        
        Returns the interface normals over the cells.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-0.5, 0.5, 0.5, 1.5))
           >>> v = 1 / Numeric.sqrt(2)
           >>> answer = Numeric.array((((0, 0), (0, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0), (v, v)),
           ...                         ((v, v), (0, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0), (0, 0))))
           >>> Numeric.allclose(distanceVariable.getCellInterfaceNormals(), answer)
           1
           
        """

        N = self.mesh.getNumberOfCells()
        M = self.mesh.getMaxFacesPerCell()
        dim = self.mesh.getDim()

        valueOverFaces = Numeric.resize(Numeric.repeat(self.getCellValueOverFaces(), dim), (N, M, dim))
        interfaceNormals = Numeric.take(self.getInterfaceNormals(), self.mesh.getCellFaceIDs())
        return Numeric.where(valueOverFaces < 0, 0, interfaceNormals)

    def getInterfaceNormals(self):
        """

        Returns the normals on the boundary faces only, the other are set to zero.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-0.5, 0.5, 0.5, 1.5))
           >>> v = 1 / Numeric.sqrt(2)
           >>> answer = Numeric.array(((0, 0), (0, 0),
           ...                         (v, v), (0, 0),
           ...                         (0, 0), (0, 0),
           ...                         (0, 0), (v, v), (0, 0),
           ...                         (0, 0), (0, 0), (0, 0)))
           >>> Numeric.allclose(distanceVariable.getInterfaceNormals(), answer)
           1
           
        """
        
        N = self.mesh.getNumberOfFaces()
        M = self.mesh.getDim()
        interfaceFlag = Numeric.resize(Numeric.repeat(self.getInterfaceFlag(), M),(N, M))
        return Numeric.where(interfaceFlag, self.getLevelSetNormals(), 0)

    def getInterfaceFlag(self):
        """

        Returns 1 for faces on boundary and 0 otherwise.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-0.5, 0.5, 0.5, 1.5))
           >>> answer = Numeric.array((0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
           >>> Numeric.allclose(distanceVariable.getInterfaceFlag(), answer)
           1
           
        """
        val0 = Numeric.take(self.value, self.mesh.getAdjacentCellIDs()[0])
        val1 = Numeric.take(self.value, self.mesh.getAdjacentCellIDs()[1])
        
        return Numeric.where(val1 * val0 < 0, 1, 0)

##    def getZeroCellFlag(self):
##        """

##        Returns a one on each cell that lies on the interface.

##           >>> from fipy.meshes.grid2D import Grid2D
##           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
##           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-0.5, 0.5, 0.5, 1.5))
##           >>> answer = Numeric.array((1, 1, 1, 0))
##           >>> Numeric.allclose(distanceVariable.getInterfaceFlag(), answer)
##           1

##        """

        


    def getCellValueOverFaces(self):
        """

        Returns the cells values at the faces.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-0.5, 0.5, 0.5, 1.5))
           >>> answer = Numeric.array(((-.5, -.5, -.5, -.5),
           ...                         (.5, .5, .5, .5),
           ...                         (.5, .5, .5, .5),
           ...                         (1.5, 1.5, 1.5, 1.5)))
           >>> Numeric.allclose(distanceVariable.getCellValueOverFaces(), answer)
           1

        """
        
        M = self.mesh.getMaxFacesPerCell()
        N = self.mesh.getNumberOfCells()
        return Numeric.reshape(Numeric.repeat(self.value, M), (N, M))

    def getLevelSetNormals(self):
        """

        Return the face level set normals.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-0.5, 0.5, 0.5, 1.5))
           >>> v = 1 / Numeric.sqrt(2)
           >>> answer = Numeric.array(((0, 0), (0, 0), (v, v), (v, v), (0, 0), (0, 0),
           ...                         (0, 0), (v, v), (0, 0), (0, 0), (v, v), (0, 0)))
           >>> Numeric.allclose(distanceVariable.getLevelSetNormals(), answer)
           1
        """
        
        faceGrad = self.getGrad().getArithmeticFaceValue()
        faceGradMag = Numeric.where(faceGrad.getMag() > 1e-10,
                                    faceGrad.getMag(),
                                    1e-10)
        faceGrad = Numeric.array(faceGrad)

        ## set faceGrad zero on exteriorFaces
        dim = self.mesh.getDim()
        exteriorFaces = (self.mesh.getExteriorFaceIDs() * dim)[:,Numeric.NewAxis] + Numeric.resize(Numeric.arange(dim), (len(self.mesh.getExteriorFaces()),dim))
        Numeric.put(faceGrad, exteriorFaces, Numeric.zeros(exteriorFaces.shape,'d'))
        
        return faceGrad / faceGradMag[:,Numeric.NewAxis] 

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 






        
