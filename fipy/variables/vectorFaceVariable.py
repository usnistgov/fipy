#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vectorFaceVariable.py"
 #                                    created: 12/9/03 {3:22:07 PM} 
 #                                last update: 8/5/05 {8:39:03 PM} 
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

"""
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
    ValueError: frames are not aligned
    
vector field times vector

    >>> print vfv * (2,)
    >>> print (2,) * vfv
    >>> print vfv * (2,3)

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
    >>> print vfv * CellVariable(mesh = mesh, value = (1,2,3))
    Traceback (most recent call last):
          ...
    ValueError: frames are not aligned
    
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
    ValueError: frames are not aligned
    >>> print (2,3,4) * vfv
    Traceback (most recent call last):
        ...
    ValueError: frames are not aligned
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
    >>> print fv * vfv

"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.variables.variable import Variable
from fipy.variables.faceVariable import FaceVariable

from fipy.tools import numerix

class VectorFaceVariable(Variable):
    def __init__(self,mesh,name = '',value=0., unit = None):

	array = Numeric.zeros([mesh._getNumberOfFaces(), mesh.getDim()],'d')
	
	Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)

        self.indexAsFaceVar = {}

    def __call__(self, point = None, order = 0):
        if point != None:
            return self[self.getMesh()._getNearestCellID(point)]
        else:
            return Variable.__call__(self)
            
    def _getVariableClass(self):
	return VectorFaceVariable

    def getIndexAsFaceVariable(self, index):

        if not self.indexAsFaceVar.has_key(index):

            from faceVariable import FaceVariable
        
            class ItemAsVariable(FaceVariable):
                def __init__(self, var, index):
                    FaceVariable.__init__(self, mesh = var.getMesh())
                    self.var = self._requires(var)

                def _calcValue(self):
                    self.value = Numeric.array(self.var[:, index])

            return ItemAsVariable(self, index)

        else:
            
            return self.indexAsFaceVar[index]
            
    def dot(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.dot(a,b), other, parentClass = FaceVariable)

    def getDivergence(self):
        if not hasattr(self, 'divergence'):
            from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
            self.divergence = _AddOverFacesVariable(self.dot(self.getMesh()._getOrientedAreaProjections()))
            
        return self.divergence
        
    def _getShapeFromMesh(mesh):
        return (mesh._getNumberOfFaces(), mesh.getDim())
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)

    def _getArithmeticParentClass(self, other):
        shape = VectorFaceVariable._getObjectShape(other)
        from fipy.variables.faceVariable import FaceVariable
        if isinstance(other, FaceVariable) \
        or shape == (self.getMesh().getDim(),) \
        or shape == (self.getMesh()._getNumberOfFaces(),) \
        or shape == (1,):
            return VectorFaceVariable
        else:
            return Variable._getArithmeticParentClass(self, other)

    def _getArithmeticBaseClass(self):
        return VectorFaceVariable

	
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
