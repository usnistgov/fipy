#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 4/4/05 {3:14:20 PM} 
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

This input file again solves a 1D diffusion problem as in::
    
    $ examples/diffusion/steadyState/mesh1D/input.py
    
The difference in this example is that the solution method is explicit.
The equation used is the `ExplicitDiffusionEquation`. In this case many
steps have to be taken to reach equilibrium.

A loop is required to execute the necessary time steps:

    >>> for step in range(steps):
    ...     var.updateOld()
    ...     eqn.solve(var, boundaryConditions = boundaryConditions, dt = timeStepDuration)

The result is again tested in the same way:

    >>> Lx = nx * dx
    >>> x = mesh.getCellCenters()[:,0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> import Numeric
    >>> ## print var.allclose(analyticalArray, rtol = 1e-3, atol = 1e-3)
    >>> print var.allclose(answer)
    1

"""

__docformat__ = 'restructuredtext'

import Numeric

from fipy.meshes.grid1D import Grid1D
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.variables.cellVariable import CellVariable
from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm
from fipy.terms.transientTerm import TransientTerm

dx = 1.
nx = 50
valueLeft = 0.
valueRight = 1.
timeStepDuration = 0.2
steps = 10

mesh = Grid1D(dx = dx, nx = nx)

var = CellVariable(
    name = "concentration",
    mesh = mesh,
    value = valueLeft)

eqn = TransientTerm() - ExplicitDiffusionTerm()

boundaryConditions=(FixedValue(mesh.getFacesLeft(),valueLeft),
                    FixedValue(mesh.getFacesRight(),valueRight))


answer = Numeric.array([  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        2.04800000e-07,  6.34880000e-06,  9.13408000e-05,  8.10188800e-04,
        4.96660480e-03,  2.23737856e-02,  7.69755136e-02,  2.07879578e-01,
        4.50699674e-01,  8.01663386e-01])

if __name__ == '__main__':
    for step in range(steps):
        var.updateOld()
        eqn.solve(var, boundaryConditions = boundaryConditions, dt = timeStepDuration)
        
    import fipy.viewers
    viewer = fipy.viewers.make(vars = var)
    viewer.plot()
    raw_input('finished')

