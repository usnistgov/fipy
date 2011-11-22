#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cellVolumeAverageVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
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

__all__ = []

from fipy.tools import numerix

from fipy.variables.variable import Variable
from fipy.variables.cellVariable import CellVariable

class _CellVolumeAverageVariable(Variable):
    """
    Takes a `CellVariable` and evaluates its volume average over all the
    cells.

        >>> from fipy.meshes import Grid2D
        >>> mesh = Grid2D(nx = 2, ny = 2, dx = 2., dy = 5.)
        >>> from fipy.variables.cellVariable import CellVariable 
        >>> var = CellVariable(value = (1, 2, 3 ,4), mesh = mesh)
        >>> print _CellVolumeAverageVariable(var)
        2.5
        
    """
    def __init__(self, var):
        Variable.__init__(self, unit = var.unit)
        self.var = self._requires(var)

    def _calcValue(self):
        mesh = self.var.mesh
        volumes = CellVariable(mesh=mesh, value=mesh.cellVolumes)
        return (self.var * volumes).sum() / volumes.sum()

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
