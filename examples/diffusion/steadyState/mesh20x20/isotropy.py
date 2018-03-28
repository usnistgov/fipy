#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "isotropy.py"
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

This input file solves a steady-state 1D diffusion problem as in
`./examples/diffusion/mesh1D.py`. The difference being that it uses a tensor for
the diffusion coefficient, even though the coefficient is isotropic.

>>> from fipy import Grid2D, CellVariable, DiffusionTerm, Viewer

>>> Lx = 20
>>> mesh = Grid2D(nx=20, ny=20)
>>> x, y = mesh.cellCenters

>>> valueLeft = 0.
>>> valueRight = 1.

>>> var = CellVariable(name = "solution variable",
...                    mesh = mesh,
...                    value = valueLeft)

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

>>> DiffusionTerm(coeff=(((1., 0.),
...                       (0., 1.)),)).solve(var)

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var).plot()

>>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
>>> print var.allclose(analyticalArray, atol = 0.025)
1
"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
