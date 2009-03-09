#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "electrostatics.py"
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

r"""
The Poisson equation is a particular example of the steady-state diffusion
equation. We examine a few cases in one dimension.

.. raw:: latex

   \IndexClass{Grid1D}

..

    >>> from fipy import *

    >>> nx = 200
    >>> dx = 0.01
    >>> L = nx * dx
    >>> mesh = Grid1D(dx = dx, nx = nx)

.. raw:: latex

   Given the electrostatic potential $\phi$,
   \IndexClass{CellVariable}
   
..

    >>> potential = CellVariable(mesh=mesh, name='potential', value=0.)

.. raw:: latex

   the permittivity $\epsilon$,

..

    >>> permittivity = 1
   
.. raw:: latex

   the concentration $C_j$ of the $j^\text{th}$ component with valence
   $z_j$ (we consider only a single component \( C_\text{e}^{-} \) with
   valence with $z_{\text{e}^{-}} = -1$)
   
..

    >>> electrons = CellVariable(mesh=mesh, name='e-')
    >>> electrons.valence = -1

.. raw:: latex

   and the charge density $\rho$,

..

    >>> charge = electrons * electrons.valence
    >>> charge.name = "charge"

The dimensionless Poisson equation is

.. raw:: latex

   \begin{equation*}
   \nabla\cdot\left(\epsilon\nabla\phi\right) = -\rho = -\sum_{j=1}^n z_j C_j
   \end{equation*}
   \IndexClass{ImplicitDiffusionTerm}

..

    >>> potential.equation = ImplicitDiffusionTerm(coeff = permittivity) \
    ...                      + charge == 0

Because this equation admits an infinite number of potential profiles,
we must constrain the solution by fixing the potential at one point:
    
.. raw:: latex

   \IndexClass{FixedValue}

..

    >>> bcs = (FixedValue(faces=mesh.getFacesLeft(), value=0),)

First, we obtain a uniform charge distribution by setting a uniform concentration
of electrons

.. raw:: latex

   $C_{\text{e}^{-}} = 1$.
   
..

    >>> electrons.setValue(1.)

and we solve for the electrostatic potential

    >>> potential.equation.solve(var=potential, 
    ...                          boundaryConditions=bcs)
   
This problem has the analytical solution

.. raw:: latex

   $$\psi(x) = \frac{x^2}{2} - 2x$$

..

    >>> x = mesh.getCellCenters()[0]
    >>> analytical = CellVariable(mesh=mesh, name="analytical solution", 
    ...                           value=(x**2)/2 - 2*x)

which has been satisifactorily obtained

    >>> print potential.allclose(analytical, rtol = 2e-5, atol = 2e-5)
    1

If we are running the example interactively, we view the result

.. raw:: latex

   \IndexModule{viewers}

..

    >>> if __name__ == '__main__':
    ...     viewer = Viewer(vars=(charge, potential, analytical))
    ...     viewer.plot()
    ...     raw_input("Press any key to continue...")

.. image:: examples/diffusion/electrostatics/uniform.pdf
   :scale: 50
   :align: center
    
Next, we segregate all of the electrons to right side of the domain

.. raw:: latex

   $$ C_{\text{e}^{-}} =
   \begin{cases}
       0& \text{for $x \le L/2$,} \\
       1& \text{for $x > L/2$.}
   \end{cases} $$
    
..

    >>> x = mesh.getCellCenters()[0]
    >>> electrons.setValue(0.)
    >>> electrons.setValue(1., where=x > L / 2.)

and again solve for the electrostatic potential

    >>> potential.equation.solve(var=potential, 
    ...                          boundaryConditions=bcs)
    
which now has the analytical solution

.. raw:: latex

   $$\psi(x) =
   \begin{cases}
       -x& \text{for $x \le L/2$,} \\
       \frac{(x-1)^2}{2} - x& \text{for $x > L/2$.}
   \end{cases} $$

..
    
    >>> analytical.setValue(-x)
    >>> analytical.setValue(((x-1)**2)/2 - x, where=x > L/2)

    >>> print potential.allclose(analytical, rtol = 2e-5, atol = 2e-5)
    1
    
and again view the result

    >>> if __name__ == '__main__':
    ...     viewer.plot()
    ...     raw_input("Press any key to continue...")

.. image:: examples/diffusion/electrostatics/right.pdf
   :scale: 50
   :align: center

Finally, we segregate all of the electrons to the left side of the
domain

.. raw:: latex

   $$ C_{\text{e}^{-}} =
   \begin{cases}
       1& \text{for $x \le L/2$,} \\
       0& \text{for $x > L/2$.}
   \end{cases} $$
    
..

    >>> electrons.setValue(1.)
    >>> electrons.setValue(0., where=x > L / 2.)

and again solve for the electrostatic potential

    >>> potential.equation.solve(var=potential, 
    ...                          boundaryConditions=bcs)
    
which has the analytical solution

.. raw:: latex

   $$\psi(x) =
   \begin{cases}
       \frac{x^2}{2} - x& \text{for $x \le L/2$,} \\
       -\frac{1}{2}& \text{for $x > L/2$.}
   \end{cases} $$

We again verify that the correct equilibrium is attained
    
    >>> analytical.setValue((x**2)/2 - x)
    >>> analytical.setValue(-0.5, where=x > L / 2)

    >>> print potential.allclose(analytical, rtol = 2e-5, atol = 2e-5)
    1
    
and once again view the result

    >>> if __name__ == '__main__':
    ...     viewer.plot()

.. image:: examples/diffusion/electrostatics/left.pdf
   :scale: 50
   :align: center

"""
__docformat__ = 'restructuredtext'
 
if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
	    
    raw_input("finished")

