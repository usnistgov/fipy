#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 3/7/05 {5:31:47 PM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""
   >>> eq.solve(var,
   ...          boundaryConditions = BCs,
   ...          solver = solver)

   >>> print var.allclose(mesh.getCellCenters()[:,0], rtol = 1e-5)
   1

"""
__docformat__ = 'restructuredtext'

from fipy.meshes.grid1D import Grid1D
dx = 1.0
nx = 100000
mesh = Grid1D(dx = dx, nx = nx)

import Numeric
from fipy.variables.cellVariable import CellVariable
var = CellVariable(mesh = mesh,
                   value = Numeric.arange(nx) * dx + dx / 2.)

from fipy.solvers.linearLUSolver import LinearLUSolver
from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
from fipy.terms.nthOrderDiffusionTerm import NthOrderDiffusionTerm
from fipy.terms.transientTerm import TransientTerm
     
eq = NthOrderDiffusionTerm((1.0,1.0))
BCs = (NthOrderBoundaryCondition(mesh.getFacesLeft(), 0., 0),
       NthOrderBoundaryCondition(mesh.getFacesRight(), nx * dx, 0),
       NthOrderBoundaryCondition(mesh.getFacesLeft(), 0., 2),
       NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 2))
solver = LinearLUSolver(steps = 10)

if __name__ == '__main__':
    eq.solve(var,
             boundaryConditions = BCs,
             solver = solver)
    
    from fipy.viewers import make
    viewer = make(var)
    viewer.plot()

    raw_input("finished")
