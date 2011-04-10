#!/usr/bin/env python

## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "vertexVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # #############################################################################
 ##

from fipy.variables.meshVariable import _MeshVariable

class _VertexVariable(_MeshVariable):
    @property
    def _variableClass(self):
        return _VertexVariable
        
    def _getShapeFromMesh(mesh):
        """
        Return the shape of this variable type, given a particular mesh.
        """
        return (mesh.numberOfVertices,)
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)
    
    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        if other is None:
            return _VertexVariable
            
        return _MeshVariable._getArithmeticBaseClass(self, other)

    @property
    def globalValue(self):
        return self._getGlobalValue(self.mesh._localNonOverlappingVertexIDs, 
                                    self.mesh._globalNonOverlappingVertexIDs)

    def _getArithmeticFaceValue(self):
        if not hasattr(self, 'arithmeticFaceValue'):
            from fipy.tools import numerix
            faceValue = numerix.take(self, self.mesh.faceVertexIDs, axis=1)
            self.arithmeticFaceValue = faceValue.mean(axis=1).filled()

        return self.arithmeticFaceValue

    @property 
    def _globalNumberOfElements(self): 
        return self.mesh.globalNumberOfVertices

    @property
    def _globalOverlappingIDs(self):
        return self.mesh._globalOverlappingVertexIDs

    @property
    def _localNonOverlappingIDs(self):
        return self.mesh._localNonOverlappingVertexIDs

