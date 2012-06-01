#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "coupledCellVariable.py"
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
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix
from fipy.tools.decorators import getsetDeprecated

class _CoupledCellVariable(object):
    def __init__(self, vars):
        self.vars = vars

    @property
    def shape(self):
        return (len(self.vars) * self.mesh.numberOfCells,)
        
    @property
    def rank(self):
        return self.vars[0].rank
        
    def __repr__(self):
        return "(" + ", ".join([repr(var) for var in self.vars]) + ")"
        
    @getsetDeprecated
    def getMesh(self):
        return self.mesh

    @property
    def mesh(self):
        meshes = list(set([var.mesh for var in self.vars]))
        
        if len(meshes) == 0:
            raise Exception("There are no Meshes defined")
        elif len(meshes) > 1:
            raise Exception("All Variables must be defined on the same Mesh")
        else:
            return meshes[0]
            
    def __getitem__(self, index):
        return numerix.concatenate([numerix.array(var[index]) for var in self.vars])
        
    def __setitem__(self, index, value):
        N = self.mesh.numberOfCells
        for i, var in enumerate(self.vars):
            if numerix.shape(value) == ():
                var[index] = value
            else:
                var[index] = value[i * N:(i + 1) * N]

    @getsetDeprecated
    def getValue(self):
        return self.value
    
    @getsetDeprecated
    def setValue(self, value):
        self.value = value

    def _getValue(self):
        return numerix.concatenate([numerix.array(var.value) for var in self.vars])

    def _setValue(self, value):
        self[:] = value

    value = property(_getValue, _setValue)
      
    @getsetDeprecated
    def getGlobalValue(self):
        return self.globalValue

    @property
    def globalValue(self):
        return numerix.concatenate([numerix.array(var.globalValue) for var in self.vars])

    @getsetDeprecated
    def getNumericValue(self):
        return self.numericValue

    @property
    def numericValue(self):
        return numerix.concatenate([var.numericValue for var in self.vars])
  
    @getsetDeprecated
    def getUnit(self):
        return self.unit

    @property
    def unit(self):
        from fipy.tools.dimensions import physicalField
        return physicalField._unity
        
    def __array__(self, t=None):
        """
        Attempt to convert the `_CoupledCellVariable` to a numerix `array` object

        >>> from fipy import *
        >>> mesh = Grid1D(nx=2)
        >>> v1 = CellVariable(mesh=mesh, value=[2, 3])
        >>> v2 = CellVariable(mesh=mesh, value=[4, 5])
        >>> v = _CoupledCellVariable(vars=(v1, v2))
        >>> print numerix.allequal([2,3,4,5], numerix.array(v)) # doctest: +PROCESSOR_0
        True
        >>> v[:] = (6,7,8,9)
        >>> print v1
        [6 7]
        >>> print v2
        [8 9]
        >>> v.getsctype() == numerix.NUMERIX.obj2sctype(numerix.array(1))
        True
        
        """
        return numerix.array(self.value, t)
        
    def __neg__(self):
        return _CoupledCellVariable([-var for var in self.vars])
        
    def __abs__(self):
        return _CoupledCellVariable([abs(var) for var in self.vars])
        
    def __iter__(self):
        return iter(self.value)

    def getsctype(self, default=None):
        if not hasattr(self, 'typecode'):
            self.typecode = numerix.obj2sctype(rep=self.numericValue, default=default)        
        return self.typecode

    def ravel(self):
        return self.value.ravel()
        
    def copy(self):
        return self.__class__(vars=[var.copy() for var in self.vars])

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
