#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "input1DpoissonAllCharge.py"
 #                                    created: 1/15/04 {3:45:27 PM} 
 #                                last update: 10/5/04 {4:59:47 PM} 
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
A simple 1D Poisson problem to test the `PoissonEquation` element of
ElPhF.

The dimensionless Poisson equation is

.. raw:: latex

   $$ \nabla\cdot\left(\epsilon\nabla\psi\right) = \rho = \sum_{j=1}^N z_j C_j$$

where 

.. raw:: latex

   $\psi$ is the electrostatic potential,
   $\epsilon$  is the permittivity,
   $\rho$ is the charge density,
   $C_j$ is the concentration of the $j$th component, and
   $z_j$ is the valence of the $j$th component.
   
We test a uniform distribution of electrons with charge

.. raw:: latex
  
   $\rho = -1$ by setting $z_j = -1$ and $C_{\text{e}^{-}} = 1$ 
   and we let the permittivity $\epsilon = 1$.
   
We iterate one timestep to equilibrate

    >>> it.timestep()

This problem has the analytical solution

.. raw:: latex

   $$\psi(x) = \frac{x^2}{2} - 2x$$

We verify that the correct equilibrium is attained

    >>> x = mesh.getCellCenters()[:,0]
    >>> analyticalArray = (x**2)/2 - 2*x

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

fields['interstitials'][0].setValue(1.)

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

