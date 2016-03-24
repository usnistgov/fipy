#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "input2D.py"
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
 #
 # ###################################################################
 ##

r"""
Solves the Cahn-Hilliard problem in a 3D cube

>>> from fipy import CellVariable, Grid3D, Viewer, GaussianNoiseVariable, TransientTerm, DiffusionTerm, DefaultSolver
>>> from fipy.tools import numerix

The only difference from :mod:`examples.cahnHilliard.mesh2D` is the
declaration of ``mesh``.

>>> if __name__ == "__main__":
...     nx = ny = nz = 100
... else:
...     nx = ny = nz = 10
>>> mesh = Grid3D(nx=nx, ny=ny, nz=nz, dx=0.25, dy=0.25, dz=0.25)
>>> phi = CellVariable(name=r"$\phi$", mesh=mesh)

We start the problem with random fluctuations about :math:`\phi = 1/2`

>>> phi.setValue(GaussianNoiseVariable(mesh=mesh,
...                                    mean=0.5,
...                                    variance=0.01))

:term:`FiPy` doesn't plot or output anything unless you tell it to:

>>> if __name__ == "__main__":
...     viewer = Viewer(vars=(phi,), datamin=0., datamax=1.)

For :term:`FiPy`, we need to perform the partial derivative
:math:`\partial f/\partial \phi`
manually and then put the equation in the canonical
form by decomposing the spatial derivatives
so that each :class:`~fipy.terms.term.Term` is of a single, even order:

.. math::

   \frac{\partial \phi}{\partial t}
    = \nabla\cdot D a^2 \left[ 1 - 6 \phi \left(1 - \phi\right)\right] \nabla \phi- \nabla\cdot D \nabla \epsilon^2 \nabla^2 \phi.

:term:`FiPy` would automatically interpolate
``D * a**2 * (1 - 6 * phi * (1 - phi))``
onto the :class:`~fipy.meshes.face.Face`\s, where the diffusive flux is calculated, but we obtain
somewhat more accurate results by performing a linear interpolation from
``phi`` at :class:`~fipy.meshes.cell.Cell` centers to ``PHI`` at :class:`~fipy.meshes.face.Face` centers.
Some problems benefit from non-linear interpolations, such as harmonic or
geometric means, and :term:`FiPy` makes it easy to obtain these, too.

>>> PHI = phi.arithmeticFaceValue
>>> D = a = epsilon = 1.
>>> eq = (TransientTerm()
...       == DiffusionTerm(coeff=D * a**2 * (1 - 6 * PHI * (1 - PHI)))
...       - DiffusionTerm(coeff=(D, epsilon**2)))

Because the evolution of a spinodal microstructure slows with time, we
use exponentially increasing time steps to keep the simulation
"interesting". The :term:`FiPy` user always has direct control over the
evolution of their problem.

>>> dexp = -5
>>> elapsed = 0.
>>> if __name__ == "__main__":
...     duration = 1000.
... else:
...     duration = 1e-2
>>> while elapsed < duration:
...     dt = min(100, numerix.exp(dexp))
...     elapsed += dt
...     dexp += 0.01
...     eq.solve(phi, dt=dt, solver=DefaultSolver(precon=None))
...     if __name__ == "__main__":
...         viewer.plot()

.. image:: mesh3D.*
   :width: 90%
   :align: center
   :alt: snapshot of Cahn-Hilliard phase separation in 3D with cutaway
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
