#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "input1DpoissonLeftCharge.py"
 #                                    created: 1/15/04 {3:45:27 PM} 
 #                                last update: 10/6/04 {4:10:40 PM} 
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
The same idea as::
    
    $ examples/elphf/input1DpoissonAllCharge.py
    
but now with charge only on the left side of the domain

.. raw:: latex

   $$ C_{\text{e}^{-}} =
   \begin{cases}
       1& \text{for $x \le L/2$,} \\
       0& \text{for $x > L/2$.}
   \end{cases} $$
   
We iterate one timestep to equilibrate

    >>> it.timestep()

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
"""
__docformat__ = 'restructuredtext'
 
from fipy.meshes.grid2D import Grid2D
from fipy.iterators.iterator import Iterator
from fipy.viewers.gist1DViewer import Gist1DViewer

import fipy.models.elphf.elphf as elphf

nx = 200
dx = 0.01
L = nx * dx

parameters = {
    'potential': {
	'name': "psi",
	'permittivity': 1.
    },
    'interstitials': (
	{
	    'name': "e-",
	    'valence': -1
	},
    )
}

mesh = Grid2D(
    dx = dx,
    dy = dx,
    nx = nx,
    ny = 1)

fields = elphf.makeFields(mesh = mesh, parameters = parameters)

setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2.)
fields['interstitials'][0].setValue(1.)
fields['interstitials'][0].setValue(0.,setCells)

equations = elphf.makeEquations(
    mesh = mesh, 
    fields = fields, 
    parameters = parameters
)

it = Iterator(equations = equations)

if __name__ == '__main__':
    viewer = Gist1DViewer(vars = (fields['charge'], fields['potential']))

    viewer.plot()
	
    raw_input("press <return> to start...")

    it.timestep()
    
    viewer.plot()
	    
    raw_input("finished")

