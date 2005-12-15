#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vectorCellVariable.py"
 #                                    created: 12/9/03 {3:22:07 PM} 
 #                                last update: 10/24/05 {5:14:47 PM} 
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

__docformat__ = 'restructuredtext'

import Numeric

from fipy.variables.variable import Variable
from fipy.variables.cellVariable import CellVariable

from fipy.tools import numerix

class VectorCellVariable(Variable):
    def __init__(self,mesh,name = '',value=0., unit = None):
	array = Numeric.zeros([mesh.getNumberOfCells(),mesh.getDim()],'d')
# 	array[:] = value	
	Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)
        self.arithmeticFaceValue = None

    def __call__(self, point = None, order = 0):
        if point != None:
            return self[self.getMesh()._getNearestCellID(point)]
        else:
            return Variable.__call__(self)

    def getArithmeticFaceValue(self):
        """

        Return a `VectorFaceVariable` with values determined by the
        arithmetic mean from the neighboring cells.
        
        >>> from fipy.meshes.grid2D import Grid2D
        >>> mesh = Grid2D(dx = 1., nx = 2)
        >>> var = VectorCellVariable(mesh, value = Numeric.array(((0,0),(1,1))))
        >>> answer = Numeric.array(((0, 0), (1, 1), (0, 0), (1, 1), 
        ...                         (0, 0), (.5, .5), (1, 1)))
        >>> Numeric.allclose(answer, Numeric.array(var.getArithmeticFaceValue()))
        1
        """
        
	if self.arithmeticFaceValue is None:
	    from vectorArithmeticCellToFaceVariable import _VectorArithmeticCellToFaceVariable
	    self.arithmeticFaceValue = _VectorArithmeticCellToFaceVariable(self)

	return self.arithmeticFaceValue
	
    def _getVariableClass(self):
	return VectorCellVariable

    def dot(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: numerix.dot(a,b), other, baseClass = CellVariable)

    def getDivergence(self):
        if not hasattr(self, 'divergence'):
            from fipy.variables.vectorFaceDifferenceVariable import _VectorFaceDifferenceVariable
            self.divergence = _VectorFaceDifferenceVariable(self).getMag()
            
        return self.divergence

    def _getShapeFromMesh(mesh):
        """
        Return the shape of this variable type, given a particular mesh.
        """
        return (mesh.getNumberOfCells(), mesh.getDim())
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)

    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            return VectorCellVariable
            
        # a VectorCellVariable operating with a CellVariable, a vector, a list of scalars, or a scalar
        # will produce a VectorCellVariable
        shape = numerix.getShape(other)
        from fipy.variables.cellVariable import CellVariable
        if isinstance(other, CellVariable) \
        or shape == (self.getMesh().getDim(),) \
        or shape == (self.getMesh().getNumberOfCells(),) \
        or shape == (1,):
            return VectorCellVariable
        else:
            return Variable._getArithmeticBaseClass(self, other)

    def _verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass, rotateVariable):
        if isinstance(var1, CellVariable) and self.getMesh() == var1.getMesh():
            return self._rotateShape(op, var0, var1, var0Array, var1Array, opShape)
        else:
            return Variable._verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass, rotateVariable)

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
