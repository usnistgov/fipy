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
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
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
