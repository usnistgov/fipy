#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "input1DpoissonLeftCharge.py"
 #                                    created: 1/15/04 {3:45:27 PM} 
 #                                last update: 12/10/04 {1:56:24 PM} 
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
The same idea as `examples/elphf/input1DpoissonAllCharge.py`, again on a 1D mesh

    >>> nx = 200
    >>> dx = 0.01
    >>> L = nx * dx
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx = dx, nx = nx)
    
but now with charge only on the left side of the domain.

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

We again let the ElPhF module construct the appropriate fields

    >>> import fipy.models.elphf.elphf as elphf
    >>> fields = elphf.makeFields(mesh = mesh, 
    ...                           parameters = parameters)

We segregate the charge to the left side of the domain by setting the
concentration of electrons to

.. raw:: latex

   $$ C_{\text{e}^{-}} =
   \begin{cases}
       1& \text{for $x \le L/2$,} \\
       0& \text{for $x > L/2$.}
   \end{cases} $$
   
..

    >>> setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2.)
    >>> fields['interstitials'][0].setValue(1.)
    >>> fields['interstitials'][0].setValue(0.,setCells)

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
       \frac{x^2}{2} - x& \text{for $x \le L/2$,} \\
       -\frac{1}{2}& \text{for $x > L/2$.}
   \end{cases} $$

We verify that the correct equilibrium is attained

    >>> import Numeric
    
    >>> x = mesh.getCellCenters()[:,0]
    >>> analyticalArray = Numeric.where(x < 1, (x**2)/2 - x, -0.5)

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

