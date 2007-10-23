## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "vertexVariable.py"
 #                                     created: 5/2/07 {10:05:56 PM}
 #                                 last update: 10/22/07 {4:22:46 PM}
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
    def _getVariableClass(self):
        return _VertexVariable
        
    def _getShapeFromMesh(mesh):
        """
        Return the shape of this variable type, given a particular mesh.
        """
        return (mesh._getNumberOfVertices(),)
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)
    
    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        if other is None:
            return _VertexVariable
            
        return _MeshVariable._getArithmeticBaseClass(self, other)


