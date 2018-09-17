#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "input.py"
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

This example is a 1D steady state diffusion test case as in
`./examples/diffusion/variable/mesh2x1/input.py` with then
number of cells set to `nx = 10`.

A simple analytical answer can be used to test the result:
   >>> DiffusionTerm(coeff = diffCoeff).solve(var)
   >>> if __name__ == "__main__":
   ...     viewer = Viewer(vars = var)
   ...     viewer.plot()
   >>> x = mesh.cellCenters[0]
   >>> values = numerix.where(x < 3. * L / 4., 10 * x - 9. * L / 4., x + 18. * L / 4.)
   >>> values = numerix.where(x < L / 4., x, values)
   >>> print var.allclose(values, atol = 1e-8, rtol = 1e-8)
   1

"""

from fipy import FaceVariable, Tri2D, CellVariable, DiffusionTerm, Viewer
from fipy.tools import numerix

nx = 10
ny = 1

valueLeft = 0.
fluxRight = 1.
timeStepDuration = 1.

L = 10.

dx = L / nx
dy = 1.

mesh = Tri2D(dx, dy, nx, ny)

var = CellVariable(
    name = "solution variable",
    mesh = mesh,
    value = valueLeft)

diffCoeff = FaceVariable(mesh = mesh, value = 1.0)

x = mesh.faceCenters[0]
diffCoeff.setValue(0.1, where=(L/4. <= x) & (x < 3. * L / 4.))

var.faceGrad.constrain([[1.], [0.]], mesh.facesRight)

var.constrain(valueLeft, mesh.facesLeft)

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    raw_input('finished')
