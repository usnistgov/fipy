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

from fipy.tools import numerix

class _CoupledCellVariable:
    def __init__(self, vars):
        self.vars = vars
        
    def __repr__(self):
        return "(" + ", ".join([repr(var) for var in self.vars]) + ")"
        
    def getMesh(self):
        meshes = list(set([var.getMesh() for var in self.vars]))
        
        if len(meshes) == 0:
            raise Exception("There are no Meshes defined")
        elif len(meshes) > 1:
            raise Exception("All Variables must be defined on the same Mesh")
        else:
            return meshes[0]
            
    def __getitem__(self, index):
        return numerix.concatenate([numerix.array(var[index]) for var in self.vars])
        
    def __setitem__(self, index, value):
        N = self.getMesh().getNumberOfCells()
        for i, var in enumerate(self.vars):
            if numerix.shape(value) == ():
                var[index] = value
            else:
                var[index] = value[i * N:(i + 1) * N]

    def getValue(self):
        return numerix.concatenate([numerix.array(var.getValue()) for var in self.vars])

    def getNumericValue(self):
        return numerix.concatenate([var.getNumericValue() for var in self.vars])

    def setValue(self, value):
        self[:] = value
        
    def getUnit(self):
        from fipy.tools.dimensions import physicalField
        return physicalField._unity
        
    def __array__(self, t=None):
        """
        Attempt to convert the `_CoupledCellVariable` to a numerix `array` object
    
        >>> mesh = Grid1D(nx=2)
        >>> v1 = CellVariable(mesh=mesh, value=[2, 3])
        >>> v2 = CellVariable(mesh=mesh, value=[4, 5])
        >>> v = _CoupledCellVariable(vars=(v1, v2))
        >>> print numerix.array(v)
        [2 3 4 5]
        """
        return numerix.array(self.getValue(), t)
        
    def __neg__(self):
        return _CoupledCellVariable([-var for var in self.vars])
        
    def __abs__(self):
        return _CoupledCellVariable([abs(var) for var in self.vars])
        
    def __iter__(self):
        return iter(self.getValue())
