#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh20x20.py"
 #                                    created: 4/6/06 {10:50:18 AM}
 #                                last update: 6/2/06 {12:38:04 PM} 
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
 #  2006- 4- 6 JEG 1.0 original
 # ###################################################################
 ##

r"""
This example again solves a 1D diffusion problem as in
``examples/diffusion/mesh1D.py``, but now on a two-dimensional mesh:

.. raw:: latex

   \IndexClass{Grid2D}

..
    
    >>> nx = 20
    >>> ny = nx
    >>> dx = 1.
    >>> dy = dx
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

We create a `CellVariable` and initialize it to zero:
    
.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phi = CellVariable(name = "solution variable",
    ...                    mesh = mesh,
    ...                    value = 0.)

and then create a diffusion equation.  This is solved by default with an
iterative conjugate gradient solver.  

.. raw:: latex

   \IndexClass{TransientTerm}
   \IndexClass{ImplicitDiffusionTerm}

..

    >>> D = 1.
    >>> from fipy.terms.transientTerm import TransientTerm
    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> eq = TransientTerm() == ImplicitDiffusionTerm(coeff=D)

We apply Dirichlet boundary conditions

    >>> valueLeft = 1
    >>> valueRight = 0

to the left and right.  Neumann boundary conditions are
automatically applied to the top and bottom.

.. raw:: latex

   \IndexClass{FixedValue}

..

    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> BCs = (FixedValue(faces=mesh.getFacesLeft(), 
    ...                   value=valueLeft),
    ...        FixedValue(faces=mesh.getFacesRight(), 
    ...                   value=valueRight))
    
We create a viewer to see the results

.. raw:: latex

   \IndexModule{viewers}

..

    >>> if __name__ == '__main__':
    ...     from fipy import viewers
    ...     viewer = viewers.make(vars=phi,
    ...                           limits={'datamin': 0., 'datamax': 1.})
    ...     viewer.plot()

and solve the equation by repeatedly looping in time:

    >>> timeStepDuration = 10 * 0.9 * dx**2 / (2 * D)
    >>> steps = 10
    >>> for step in range(steps):
    ...     eq.solve(var=phi,
    ...              boundaryConditions=BCs,
    ...              dt=timeStepDuration)
    ...     if __name__ == '__main__':
    ...         viewer.plot()

.. image:: examples/diffusion/mesh20x20transient.pdf
   :scale: 50
   :align: center

We can again test against the analytical solution 

.. raw:: latex

   $\phi = 1 - \erf(x/2\sqrt{D t})$.
   \IndexModule{numerix}
   \IndexSoftware{SciPy}
   \IndexFunction{sqrt}

..

    >>> x = mesh.getCellCenters()[0]
    >>> t = timeStepDuration * steps
    >>> from fipy.tools.numerix import sqrt

    >>> phiAnalytical = CellVariable(name="analytical value",
    ...                              mesh=mesh)

    >>> try:
    ...     from scipy.special import erf
    ...     phiAnalytical.setValue(1 - erf(x / (2 * sqrt(D * t))))
    ... except ImportError:
    ...     print "The SciPy library is not available to test the solution to \
    ... the transient diffusion equation"

    >>> print phi.allclose(phiAnalytical, atol = 4e-2)
    1

    >>> if __name__ == '__main__':
    ...     raw_input("Implicit transient diffusion. Press <return> to proceed...")

-----

We can also solve the steady-state problem directly

    >>> ImplicitDiffusionTerm().solve(var=phi, 
    ...                               boundaryConditions = BCs)
    >>> if __name__ == '__main__':
    ...     viewer.plot()

.. image:: examples/diffusion/mesh20x20steadyState.pdf
   :scale: 50
   :align: center

and test the result against the expected linear composition profile:

    >>> L = nx * dx
    >>> x = mesh.getCellCenters()[0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / L
    >>> print phi.allclose(analyticalArray, rtol = 2e-10, atol = 2e-10)
    1
    
    >>> if __name__ == '__main__':
    ...     raw_input("Implicit steady-state diffusion. Press <return> to proceed...")

"""

__docformat__ = 'restructuredtext'

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

