#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mixedelement.py"
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

This input file again solves an explicit 1D diffusion problem as in
`./examples/diffusion/mesh1D.py` but on a mesh with both square and triangular
elements. The term used is the `ExplicitDiffusionTerm`. In this case many steps
have to be taken to reach equilibrium. The `timeStepDuration` parameter specifies
the size of each time step and `steps` is the number of time steps.

>>> from fipy import CellVariable, Grid2D, Tri2D, TransientTerm, ExplicitDiffusionTerm, Viewer
>>> from fipy.tools import numerix

>>> dx = 1.
>>> dy = 1.
>>> nx = 10

>>> valueLeft = 0.
>>> valueRight = 1.
>>> D = 1.
>>> timeStepDuration = 0.01 * dx**2 / (2 * D)
>>> if __name__ == '__main__':
...     steps = 10000
... else:
...     steps = 10

>>> gridMesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=2)
>>> triMesh = Tri2D(dx=dx, dy=dy, nx=nx, ny=1) + ((dx*nx,), (0,))
>>> otherGridMesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=1) + ((dx*nx,), (1,))
>>> bigMesh = gridMesh + triMesh + otherGridMesh

>>> L = dx * nx * 2

>>> var = CellVariable(name="concentration",
...                    mesh=bigMesh,
...                    value=valueLeft)

>>> eqn = TransientTerm() == ExplicitDiffusionTerm(coeff=D)

>>> var.constrain(valueLeft, where=bigMesh.facesLeft)
>>> var.constrain(valueRight, where=bigMesh.facesRight)

In a semi-infinite domain, the analytical solution for this transient diffusion
problem is given by :math:`\phi = 1 - \erf((L - x)/2\sqrt{D t})`, which is a
reasonable approximation at early times. At late times, the solution is just a
straight line. If the :term:`SciPy` library is available, the result is tested
against the expected profile:

>>> x = bigMesh.cellCenters[0]
>>> t = timeStepDuration * steps

>>> if __name__ == '__main__':
...     varAnalytical = valueLeft + (valueRight - valueLeft) * x / L
...     atol = 0.2
... else:
...     try:
...         from scipy.special import erf # doctest: +SCIPY
...         varAnalytical = 1 - erf((L - x) / (2 * numerix.sqrt(D * t))) # doctest: +SCIPY
...         atol = 0.03
...     except ImportError:
...         print "The SciPy library is not available to test the solution to \
...           the transient diffusion equation"

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var)

>>> for step in range(steps):
...     eqn.solve(var, dt=timeStepDuration)
...     if (step % 100) == 0 and (__name__ == '__main__'):
...         viewer.plot()

We check the answer against the analytical result

>>> print var.allclose(varAnalytical, atol=atol) # doctest: +SCIPY
1

>>> if __name__ == '__main__':
...     viewer.plot()
...     raw_input('finished')
"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
