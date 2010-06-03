#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "faceVariable.py"
 #
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from fipy.variables.meshVariable import _MeshVariable
from fipy.tools import numerix

class FaceVariable(_MeshVariable):
    def _getVariableClass(self):
        return FaceVariable

    def _getShapeFromMesh(mesh):
        """
        Return the shape of this variable type, given a particular mesh.
        """
        return (mesh._getNumberOfFaces(),)
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)
        
    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        if other is None:
            return FaceVariable
            
        return _MeshVariable._getArithmeticBaseClass(self, other)

    def copy(self):
        return self._getArithmeticBaseClass()(mesh=self.mesh, 
                                              name=self.name + "_copy", 
                                              value=self.getValue())

    def getGlobalValue(self):
        return self._getGlobalValue(self.mesh._getLocalNonOverlappingFaceIDs(), 
                                    self.mesh._getGlobalNonOverlappingFaceIDs())

    def setValue(self, value, unit = None, where = None):
        _MeshVariable.setValue(self, value=self._globalToLocalValue(value), unit=unit, where=where)

    def getDivergence(self):
        """
            >>> from fipy.meshes.grid2D import Grid2D
            >>> from fipy.variables.cellVariable import CellVariable
            >>> mesh = Grid2D(nx=3, ny=2)
            >>> var = CellVariable(mesh=mesh, value=range(3*2))
            >>> print var.getFaceGrad().getDivergence()
            [ 4.  3.  2. -2. -3. -4.]
            
        """
        if not hasattr(self, 'divergence'):
            from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
            self.divergence = _AddOverFacesVariable(self.dot(self.getMesh()._getOrientedAreaProjections()))
            
        return self.divergence

    def _getGlobalNumberOfElements(self):
        return self.mesh.globalNumberOfFaces
        
    def _getGlobalOverlappingIDs(self):
        return self.mesh._getGlobalOverlappingFaceIDs()
        
    def _getLocalNonOverlappingIDs(self):
        return self.mesh._getLocalNonOverlappingFaceIDs()

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

