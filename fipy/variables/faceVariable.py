#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "faceVariable.py"
 #                                    created: 12/9/03 {1:58:27 PM} 
 #                                last update: 1/17/06 {11:17:03 AM} 
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

import Numeric

from fipy.variables.variable import Variable
from fipy.tools import numerix

class FaceVariable(Variable):
    def __init__(self, mesh, name = '', value=0., unit = None):
        if value is None:
            array = None
        else:
            array = Numeric.zeros(self._getShapeFromMesh(mesh),'d')
# 	array[:] = value
	Variable.__init__(self,mesh = mesh, name = name, value = value, unit = unit, array = array)

    def setValue(self, value, faces = (), unit = None, where = None):
        if faces == ():
            Variable.setValue(self, value, unit = unit, where = where)
        else:
            import warnings
            warnings.warn("'where' should be used instead of 'faces'", DeprecationWarning, stacklevel=2)
            for face in faces:
                self[face.getID()] = value

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
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            return FaceVariable
            
        # a FaceVariable operating with a vector will produce a VectorFaceVariable
        if numerix.getShape(other) == (self.getMesh().getDim(),):
            from fipy.variables.vectorFaceVariable import VectorFaceVariable
            return VectorFaceVariable
        else:
            return Variable._getArithmeticBaseClass(self, other)

    def _verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass):
        from fipy.variables.vectorFaceVariable import VectorFaceVariable
        if isinstance(var1, VectorFaceVariable) and self.getMesh() == var1.getMesh():
            return self._rotateShape(op, var1, var0, var1Array, var0Array, opShape)
        else:
            return Variable._verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass)
