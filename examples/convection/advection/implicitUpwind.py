#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "impicitUpwind.py"
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

"""
This example shows the failure of advecting a square pulse with a first
order implicit upwind scheme.
"""

from fipy import CellVariable, Grid1D, TransientTerm, PowerLawConvectionTerm, LinearLUSolver, Viewer
from fipy.tools import numerix

valueLeft = 0.
valueRight = 0.
L = 10.
nx = 400
dx = L / nx
cfl = 0.01
velocity = 1.
timeStepDuration = cfl * dx / abs(velocity)
steps = 1000

mesh = Grid1D(dx = dx, nx = nx)

startingArray = numerix.zeros(nx, 'd')
startingArray[50:90] = 1.

var = CellVariable(
    name = "advection variable",
    mesh = mesh,
    value = startingArray)

var.constrain(valueLeft, mesh.facesLeft)
var.constrain(valueRight, mesh.facesRight)

eq = TransientTerm() - PowerLawConvectionTerm(coeff = (velocity,))

if __name__ == '__main__':

    viewer = Viewer(vars=(var,))
    viewer.plot()
    raw_input("press key to continue")
    for step in range(steps):
        eq.solve(var,
                 dt = timeStepDuration,
                 solver = LinearLUSolver(tolerance = 1.e-15))
        viewer.plot()
    viewer.plot()
    raw_input('finished')
