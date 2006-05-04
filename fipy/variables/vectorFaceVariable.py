#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vectorFaceVariable.py"
 #                                    created: 12/9/03 {3:22:07 PM} 
 #                                last update: 5/4/06 {7:58:11 AM} 
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
from fipy.variables.faceVariable import FaceVariable

from fipy.tools import numerix

class VectorFaceVariable(Variable):
    def __init__(self,mesh,name = '',value=0., unit = None):
        if value is None:
            array = None
        else:
            array = Numeric.zeros(self._getShapeFromMesh(mesh),'d')
	
	Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)

    def __call__(self, point = None, order = 0):
        if point != None:
            return self[self.getMesh()._getNearestCellID(point)]
        else:
            return Variable.__call__(self)
            
    def _getVariableClass(self):
	return VectorFaceVariable

    def dot(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.dot(a,b), other, baseClass = FaceVariable)

    def getDivergence(self):
        if not hasattr(self, 'divergence'):
            from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
            self.divergence = _AddOverFacesVariable(self.dot(self.getMesh()._getOrientedAreaProjections()))
            
        return self.divergence
        
    def _getShapeFromMesh(mesh):
        """
        Return the shape of this variable type, given a particular mesh.
        """
        return (mesh._getNumberOfFaces(), mesh.getDim())
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)

    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            return VectorFaceVariable
            
        # a VectorFaceVariable operating with a FaceVariable, a vector, a list of scalars, or a scalar
        # will produce a VectorFaceVariable
        shape = numerix.getShape(other)
        from fipy.variables.faceVariable import FaceVariable
        if isinstance(other, FaceVariable) \
        or shape == (self.getMesh().getDim(),) \
        or shape == (self.getMesh()._getNumberOfFaces(),) \
        or shape == (1,):
            return VectorFaceVariable
        else:
            return Variable._getArithmeticBaseClass(self, other)

    def _verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass, rotateShape = True):
        if isinstance(var1, FaceVariable) and self.getMesh() == var1.getMesh():
            return self._rotateShape(op, var0, var1, var0Array, var1Array, opShape)
        else:
            return Variable._verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass, rotateShape)

    def _testArithmetic(self):
        r"""
            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(nx = 3)
            >>> vfv = VectorFaceVariable(mesh = mesh, value = ((0,),(1,),(2,),(3,)))
            
        vector field times scalar

            >>> print vfv * 3
            [[ 0.,]
             [ 3.,]
             [ 6.,]
             [ 9.,]]
            >>> print 3 * vfv
            [[ 0.,]
             [ 3.,]
             [ 6.,]
             [ 9.,]]

        vector field times scalar variable

            >>> from fipy.variables.variable import Variable
            >>> print vfv * Variable(value = 3)
            [[ 0.,]
             [ 3.,]
             [ 6.,]
             [ 9.,]]
            >>> print Variable(value = 3) * vfv
            [[ 0.,]
             [ 3.,]
             [ 6.,]
             [ 9.,]]

        vector field times vector field

            >>> print vfv * vfv
            [[ 0.,]
             [ 1.,]
             [ 4.,]
             [ 9.,]]

        vector field times cell centered field

            >>> from fipy.variables.cellVariable import CellVariable
            >>> print vfv * CellVariable(mesh = mesh, value = (1,2,3))
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            
        .. note:: 
            
           Older error message was:: 
            
               TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            
        vector field times vector

            >>> print vfv * (2,)
            [[ 0.,]
             [ 2.,]
             [ 4.,]
             [ 6.,]]
            >>> print (2,) * vfv
            [[ 0.,]
             [ 2.,]
             [ 4.,]
             [ 6.,]]
            >>> print vfv * (2,3)
            Traceback (most recent call last):
                ...
            TypeError: can't multiply sequence to non-int

        How about 2D meshes?

            >>> from fipy.meshes.grid2D import Grid2D
            >>> mesh = Grid2D(nx = 3, ny = 1)
            >>> vfv = VectorFaceVariable(mesh = mesh, value = ((0,1),(1,2),(2,3),(3,4),(1,3),(2,4),(3,5),(6,9),(2,6),(1,3)))
            
        vector field times scalar

            >>> print vfv * 3
            [[  0.,  3.,]
             [  3.,  6.,]
             [  6.,  9.,]
             [  9., 12.,]
             [  3.,  9.,]
             [  6., 12.,]
             [  9., 15.,]
             [ 18., 27.,]
             [  6., 18.,]
             [  3.,  9.,]]
            >>> print 3 * vfv
            [[  0.,  3.,]
             [  3.,  6.,]
             [  6.,  9.,]
             [  9., 12.,]
             [  3.,  9.,]
             [  6., 12.,]
             [  9., 15.,]
             [ 18., 27.,]
             [  6., 18.,]
             [  3.,  9.,]]
             
        vector field times scalar variable

            >>> from fipy.variables.variable import Variable
            >>> print vfv * Variable(value = 3)
            [[  0.,  3.,]
             [  3.,  6.,]
             [  6.,  9.,]
             [  9., 12.,]
             [  3.,  9.,]
             [  6., 12.,]
             [  9., 15.,]
             [ 18., 27.,]
             [  6., 18.,]
             [  3.,  9.,]]
            >>> print Variable(value = 3) * vfv
            [[  0.,  3.,]
             [  3.,  6.,]
             [  6.,  9.,]
             [  9., 12.,]
             [  3.,  9.,]
             [  6., 12.,]
             [  9., 15.,]
             [ 18., 27.,]
             [  6., 18.,]
             [  3.,  9.,]]

        vector field times vector field

            >>> print vfv * vfv
            [[  0.,  1.,]
             [  1.,  4.,]
             [  4.,  9.,]
             [  9., 16.,]
             [  1.,  9.,]
             [  4., 16.,]
             [  9., 25.,]
             [ 36., 81.,]
             [  4., 36.,]
             [  1.,  9.,]]

        vector field times cell centered field

            >>> from fipy.variables.cellVariable import CellVariable
            >>> cv = CellVariable(mesh = mesh, value = (1,2,3))
            >>> print vfv * cv
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            
        .. note:: 
            
           Older error message was::
            
               TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            
        vector field times vector

            >>> print vfv * (2,3)
            [[  0.,  3.,]
             [  2.,  6.,]
             [  4.,  9.,]
             [  6., 12.,]
             [  2.,  9.,]
             [  4., 12.,]
             [  6., 15.,]
             [ 12., 27.,]
             [  4., 18.,]
             [  2.,  9.,]]
            >>> print (2,3) * vfv
            [[  0.,  3.,]
             [  2.,  6.,]
             [  4.,  9.,]
             [  6., 12.,]
             [  2.,  9.,]
             [  4., 12.,]
             [  6., 15.,]
             [ 12., 27.,]
             [  4., 18.,]
             [  2.,  9.,]]
            >>> print vfv * (2,3,4)
            Traceback (most recent call last):
                ...
            TypeError: can't multiply sequence to non-int
            >>> print (2,3,4) * vfv
            Traceback (most recent call last):
                ...
            TypeError: can't multiply sequence to non-int
            >>> print vfv * Variable(value = (2,3))
            [[  0.,  3.,]
             [  2.,  6.,]
             [  4.,  9.,]
             [  6., 12.,]
             [  2.,  9.,]
             [  4., 12.,]
             [  6., 15.,]
             [ 12., 27.,]
             [  4., 18.,]
             [  2.,  9.,]]
            >>> print Variable(value = (2,3)) * vfv
            [[  0.,  3.,]
             [  2.,  6.,]
             [  4.,  9.,]
             [  6., 12.,]
             [  2.,  9.,]
             [  4., 12.,]
             [  6., 15.,]
             [ 12., 27.,]
             [  4., 18.,]
             [  2.,  9.,]]

        vector field times scalar field

            >>> from fipy.variables.faceVariable import FaceVariable
            >>> fv = FaceVariable(mesh = mesh, value = (0,1,2,3,4,5,6,7,8,9))
            >>> print vfv * fv
            [[  0.,  0.,]
             [  1.,  2.,]
             [  4.,  6.,]
             [  9., 12.,]
             [  4., 12.,]
             [ 10., 20.,]
             [ 18., 30.,]
             [ 42., 63.,]
             [ 16., 48.,]
             [  9., 27.,]]
            >>> print fv * vfv
            [[  0.,  0.,]
             [  1.,  2.,]
             [  4.,  6.,]
             [  9., 12.,]
             [  4., 12.,]
             [ 10., 20.,]
             [ 18., 30.,]
             [ 42., 63.,]
             [ 16., 48.,]
             [  9., 27.,]]

        """
        pass

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
