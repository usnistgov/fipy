#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "input1DpoissonRightCharge.py"
 #                                    created: 1/15/04 {3:45:27 PM} 
 #                                last update: 12/10/04 {1:52:37 PM} 
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
 #  2004-01-15 JEG 1.0 original
 # ###################################################################
 ##

r"""
A simple problem to test the `PoissonEquation` element of
ElPhF on a 1D mesh

    >>> nx = 200
    >>> dx = 0.01
    >>> L = nx * dx
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx = dx, nx = nx)

The dimensionless Poisson equation is

.. raw:: latex

   $$ \nabla\cdot\left(\epsilon\nabla\psi\right) = -\rho = -\sum_{j=1}^n z_j C_j$$

where 

.. raw:: latex

   $\psi$ is the electrostatic potential,
   $\epsilon$  is the permittivity,
   $\rho$ is the charge density,
   $C_j$ is the concentration of the $j^\text{th}$ component, and
   $z_j$ is the valence of the $j^\text{th}$ component.
   
..

We examine a fixed distribution of electrons

.. raw:: latex
  
   with $z_{\text{e}^{-}} = -1$ and and we let the permittivity $\epsilon = 1$.
   
In the ElPhF construction, electrons are treated as interstitial elements, 
which can diffuse freely without displacing other components

    >>> parameters = {
    ...     'potential': {
    ...         'name': "psi",
    ...         'permittivity': 1.,
    ...     },
    ...     'interstitials': (
    ...         {
    ...             'name': "e-",
    ...             'valence': -1,
    ...             'diffusivity': 0
    ...         },
    ...     )
    ... }

We have set the diffusivity of electrons to zero to keep them from moving due 
to electromigration.

We again let the ElPhF module construct the appropriate fields

    >>> import fipy.models.elphf.elphf as elphf
    >>> fields = elphf.makeFields(mesh = mesh, 
    ...                           parameters = parameters)

We segregate all of the electrons to one side of the domain

.. raw:: latex

   $$ C_{\text{e}^{-}} =
   \begin{cases}
       0& \text{for $x \le L/2$,} \\
       1& \text{for $x > L/2$.}
   \end{cases} $$
    
..

    >>> setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2.)
    >>> fields['interstitials'][0].setValue(0.)
    >>> fields['interstitials'][0].setValue(1.,setCells)

and iterate one implicit timestep to equilibrate the electrostatic potential

    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> bcs = (FixedValue(faces = mesh.getFacesLeft(), value = 0),)
    
    >>> from fipy.models.elphf.poissonEquation import factory
    >>> poisson = factory.make(fields, parameters['potential'])
    >>> poisson.solve(var = fields['potential'], 
    ...               boundaryConditions = bcs)
    
This problem has the analytical solution

.. raw:: latex

   $$\psi(x) =
   \begin{cases}
       -x& \text{for $x \le L/2$,} \\
       \frac{(x-1)^2}{2} - x& \text{for $x > L/2$.}
   \end{cases} $$

We verify that the correct equilibrium is attained
    
    >>> x = mesh.getCellCenters()[:,0]

    >>> import Numeric
    >>> analyticalArray = Numeric.where(x < 1, -x, ((x-1)**2)/2 - x)

    >>> fields['potential'].allclose(analyticalArray, rtol = 2e-5, atol = 2e-5)
    1
    
If we are running the example interactively, we view the result

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.gist1DViewer import Gist1DViewer
    ...     viewer = Gist1DViewer(vars = (fields['charge'], fields['potential']))
    ...     viewer.plot()
"""
__docformat__ = 'restructuredtext'
 
if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
	    
    raw_input("finished")

