#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vectorCellVariable.py"
 #                                    created: 12/9/03 {3:22:07 PM} 
 #                                last update: 10/19/04 {12:22:48 PM} 
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

import fipy.tools.array as array

class VectorCellVariable(Variable):
    def __init__(self,mesh,name = '',value=0., unit = None):
	array = Numeric.zeros([len(mesh.getCells()),mesh.getDim()],'d')
# 	array[:] = value	
	Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)
        self.arithmeticFaceValue = None

    def getArithmeticFaceValue(self):
        """

        Return a `VectorFaceVariable` with values determined by the
        arithmetic mean from the neighboring cells.
        
        >>> from fipy.meshes.grid2D import Grid2D
        >>> mesh = Grid2D(dx = 1., dy = 1, nx = 2, ny = 1)
        >>> var = VectorCellVariable(mesh, value = Numeric.array(((0,0),(1,1))))
        >>> answer = Numeric.array(((0, 0), (1, 1), (0, 0), (1, 1), (0, 0), (.5, .5), (1, 1)))
        >>> Numeric.allclose(answer, Numeric.array(var.getArithmeticFaceValue()))
        1
        """
        
	if self.arithmeticFaceValue is None:
	    from vectorArithmeticCellToFaceVariable import VectorArithmeticCellToFaceVariable
	    self.arithmeticFaceValue = VectorArithmeticCellToFaceVariable(self)

	return self.arithmeticFaceValue
	
    def getVariableClass(self):
	return VectorCellVariable

    def dot(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: array.dot(a,b), other, parentClass = CellVariable)
	

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
