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

    >>> from fipy import *

The only difference from ``examples.cahnHilliard.mesh1000x1000`` is the
declaration of `mesh`.

    >>> mesh = Grid3D(nx=100, ny=100, nz=100, dx=0.25, dy=0.25, dz=0.25)
    >>> phi = CellVariable(name=r"$\phi$", mesh=mesh)

We start the problem with random fluctuations about $\phi = 1/2$

    >>> phi.setValue(GaussianNoiseVariable(mesh=mesh,
    ...                                    mean=0.5,
    ...                                    variance=0.01))

FiPy doesn't plot or output anything unless you tell it to:

    >>> if __name__ == "__main__":
    ...     viewer = Viewer(vars=(phi,), datamin=0., datamax=1.)

For FiPy, we need to perform the partial derivative 

.. raw:: latex

   $\partial f/\partial \phi$ 
   
manually and then put the equation in the canonical
form by decomposing the spatial derivatives
so that each `Term` is of a single, even order:
    
.. raw:: latex

   $$\frac{\partial \phi}{\partial t}
    = \nabla\cdot D a^2 \left[ 1 - 6 \phi \left(1 - \phi\right)\right] \nabla \phi- \nabla\cdot D \nabla \epsilon^2 \nabla^2 \phi.$$

FiPy would automatically interpolate
`D * a**2 * (1 - 6 * phi * (1 - phi))`
onto the `Face`\s, where the diffusive flux is calculated, but we obtain
somewhat more accurate results by performing a linear interpolation from
`phi` at `Cell` centers to `PHI` at `Face` centers.
Some problems benefit from non-linear interpolations, such as harmonic or
geometric means, and FiPy makes it easy to obtain these, too.

    >>> PHI = phi.getArithmeticFaceValue()
    >>> D = a = epsilon = 1.
    >>> eq = (TransientTerm()
    ...       == DiffusionTerm(coeff=D * a**2 * (1 - 6 * PHI * (1 - PHI)))
    ...       - DiffusionTerm(coeff=(D, epsilon**2)))

Because the evolution of a spinodal microstructure slows with time, we
use exponentially increasing time steps to keep the simulation
"interesting". The FiPy user always has direct control over the
evolution of their problem.

    >>> dexp = -5
    >>> elapsed = 0.
    >>> while elapsed < 1000.:
    ...     dt = min(100, exp(dexp))
    ...     elapsed += dt
    ...     dexp += 0.01
    ...     eq.solve(phi, dt=dt)
    ...     if __name__ == "__main__":
    ...         viewer.plot()
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    
    raw_input('finished')


