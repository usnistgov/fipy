#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "gmshinput.py"
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

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. In order to test
the non-orthogonality error, this uses a SkewedGrid2D, which is a
Grid2D with each interior vertex moved in a random direction.

"""

if __name__ == '__main__':
    import sys
    import os

    from fipy import SkewedGrid2D, CellVariable, Viewer, DiffusionTerm
    from fipy.tools import numerix

    valueLeft = 0.
    valueRight = 1.

    mesh = SkewedGrid2D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, rand = 0.1)

    var = CellVariable(name = "solution variable",
                       mesh = mesh,
                       value = valueLeft)

    viewer = Viewer(vars = var)

    var.constrain(valueLeft, mesh.facesLeft)
    var.constrain(valueRight, mesh.facesRight)

    DiffusionTerm().solve(var)

    varArray = numerix.array(var)
    x = mesh.cellCenters[0]
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
    errorArray = varArray - analyticalArray
    errorVar = CellVariable(name = 'absolute error',
                            mesh = mesh,
                            value = abs(errorArray))
    errorViewer = Viewer(vars = errorVar)

    NonOrthoVar = CellVariable(name = "non-orthogonality",
                               mesh = mesh,
                               value = mesh._nonOrthogonality)
    NOViewer = Viewer(vars = NonOrthoVar)
    viewer.plot()
    NOViewer.plot()

    raw_input("finished")
