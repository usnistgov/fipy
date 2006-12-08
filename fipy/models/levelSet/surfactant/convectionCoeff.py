#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "convectionCoeff.py"
 #                                    created: 7/28/04 {10:39:23 AM} 
 #                                last update: 11/16/06 {2:01:26 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools.numerix import MA
from fipy.tools import numerix

from fipy.tools import numerix
from fipy.tools import vector

from fipy.variables.vectorFaceVariable import VectorFaceVariable

class _ConvectionCoeff(VectorFaceVariable):
    """
    
    Convection coefficient for the `ConservativeSurfactantEquation`.
    The coeff only has a value for a negative `distanceVar`.

    """

    def __init__(self, distanceVar):
        """
        
        Simple one dimensional test:

        
           >>> from fipy.variables.cellVariable import CellVariable
           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(nx = 3, ny = 1, dx = 1., dy = 1.)
           >>> from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
           >>> distanceVar = DistanceVariable(mesh, value = (-.5, .5, 1.5))
           >>> answer = numerix.zeros((mesh._getNumberOfFaces(),2),'d')
           >>> answer[7,0] = -1
           >>> print numerix.allclose(_ConvectionCoeff(distanceVar), answer)
           1

        Change the dimensions:

           >>> mesh = Grid2D(nx = 3, ny = 1, dx = .5, dy = .25)
           >>> distanceVar = DistanceVariable(mesh, value = (-.25, .25, .75))
           >>> answer[7,0] = -.5
           >>> print numerix.allclose(_ConvectionCoeff(distanceVar), answer)
           1

        Two dimensional example:

           >>> mesh = Grid2D(nx = 2, ny = 2, dx = 1., dy = 1.)
           >>> distanceVar = DistanceVariable(mesh, value = (-1.5, -.5, -.5, .5))
           >>> answer = numerix.zeros((mesh._getNumberOfFaces(),2),'d')
           >>> answer[2,1] = -.5
           >>> answer[3,1] = -1
           >>> answer[7,0] = -.5
           >>> answer[10,0] = -1
           >>> print numerix.allclose(_ConvectionCoeff(distanceVar), answer)
           1

        Larger grid:

           >>> mesh = Grid2D(nx = 3, ny = 3, dx = 1., dy = 1.)
           >>> distanceVar = DistanceVariable(mesh, value = (1.5, .5 , 1.5,
           ...                                           .5 , -.5, .5 ,
           ...                                           1.5, .5 , 1.5))
           >>> answer = numerix.zeros((mesh._getNumberOfFaces(),2), 'd')
           >>> answer[4,1] = .25
           >>> answer[7,1] = -.25
           >>> answer[7,1] = -.25
           >>> answer[17,0] = .25
           >>> answer[18,0] = -.25
           >>> print numerix.allclose(_ConvectionCoeff(distanceVar), answer)
           1
           
        """
        
        VectorFaceVariable.__init__(self, distanceVar.getMesh(), name = 'surfactant convection')
        self.distanceVar = self._requires(distanceVar)

    def _calcValue(self):

        Ncells = self.mesh.getNumberOfCells()
        Nfaces = self.mesh._getNumberOfFaces()
        M = self.mesh._getMaxFacesPerCell()
        dim = self.mesh.getDim()
        cellFaceIDs = self.mesh._getCellFaceIDs()
     
        faceNormalAreas = self.distanceVar._getLevelSetNormals() * self.mesh._getFaceAreas()[:,numerix.NewAxis]

        cellFaceNormalAreas = numerix.array(numerix.MAtake(faceNormalAreas, cellFaceIDs).filled(fill_value = 0))
        norms = numerix.array(MA.array(self.mesh._getCellNormals()).filled(fill_value = 0))
        
        alpha = numerix.dot(cellFaceNormalAreas, norms, axis = 2)
        alpha = numerix.where(alpha > 0, alpha, 0)

        alphasum = numerix.sum(alpha, axis = 1)
        alphasum += (alphasum < 1e-10) * 1e-10
        alpha = alpha / alphasum[:,numerix.NewAxis]

        phi = numerix.reshape(numerix.repeat(self.distanceVar, M), (Ncells, M))
        alpha = numerix.where(phi > 0., 0, alpha)

        volumes = numerix.array(self.mesh.getCellVolumes())
        alpha = alpha[:,:,numerix.NewAxis] * volumes[:,numerix.NewAxis,numerix.NewAxis] * norms

        value = numerix.zeros(Nfaces * dim,'d')

        cellFaceIDs = (cellFaceIDs.flat * dim)[:,numerix.NewAxis] + numerix.resize(numerix.arange(dim), (len(cellFaceIDs.flat),dim))
        
        vector._putAddPy(value, cellFaceIDs.flat, alpha.flat, mask = MA.array(cellFaceIDs).flat.mask())

        value = numerix.reshape(value, (Nfaces, dim))

        return -value / self.mesh._getFaceAreas()[:,numerix.NewAxis]

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
