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

r"""

To run this example from the base fipy directory type::

    $ python examples/diffusion/steadyState/mesh1D/tri2Dinput.py

at the command line. A contour plot should appear and the word `finished`
in the terminal.

This example is similar to the example found in
:mod:`examples.diffusion.steadyState.mesh1D.input`, however, the `mesh` is a
`Tri2D` object rather than a `Grid2D` object.

Here, one time step is executed to implicitly find the steady state
solution.

    >>> DiffusionTerm().solve(var)

To test the solution, the analytical result is required. The `x`
coordinates from the mesh are gathered and the length of the domain,
`Lx`, is calculated.  An array, `analyticalArray`, is calculated to
compare with the numerical result,

    >>> x = mesh.cellCenters[0]
    >>> Lx = nx * dx
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx

Finally the analytical and numerical results are compared with a
tolerance of `1e-10`.

    >>> print var.allclose(analyticalArray)
    1

"""

__docformat__ = 'restructuredtext'

from fipy import Tri2D, CellVariable, DiffusionTerm, Viewer

nx = 50
dx = 1.

mesh = Tri2D(dx = dx, nx = nx)

valueLeft = 0.
valueRight = 1.
var = CellVariable(name = "solution-variable", mesh = mesh, value = valueLeft)

var.constrain(valueLeft, mesh.facesLeft)
var.constrain(valueRight, mesh.facesRight)

if __name__ == '__main__':
    DiffusionTerm().solve(var)
    viewer = Viewer(vars=var)
    viewer.plot()
    x = mesh.cellCenters[0]
    Lx = nx * dx
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    print var.allclose(analyticalArray)
    raw_input("finished")
