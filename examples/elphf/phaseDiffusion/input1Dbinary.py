#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/6/04 {4:45:45 PM} 
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

r"""
This example combines a 1D phase field problem, as given in::
    
    $ examples/elphf/input1Dphase.py
    
with a binary diffusion problem.

We start with a sharp phase boundary

.. raw:: latex

   $$ \xi =
   \begin{cases}
       1& \text{for $x \le L/2$,} \\
       0& \text{for $x > L/2$,}
   \end{cases} $$
   
and with a uniform concentration field

.. raw:: latex

   $C_1 = 0.5$.

This problem does not have an analytical solution, so after iterating to
equilibrium

    >>> for i in range(40):
    ...     it.timestep()

we confirm that the far-field phases have remained separated

    >>> ends = Numeric.take(fields['phase'], (0,-1))
    >>> Numeric.allclose(ends, (1.0, 0.0), rtol = 2e-3, atol = 2e-3)
    1
    
and that the concentration field has appropriately segregated into solute
rich and solute poor phases

    >>> ends = Numeric.take(fields['substitutionals'][0], (0,-1))
    >>> Numeric.allclose(ends, (0.7, 0.3), rtol = 2e-3, atol = 2e-3)
    1
"""
__docformat__ = 'restructuredtext'

import Numeric

## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.gist1DViewer import Gist1DViewer
from fipy.iterators.iterator import Iterator

import fipy.models.elphf.elphf as elphf

nx = 400
dx = 0.01
L = nx * dx

mesh = Grid2D(
    dx = dx,
    dy = dx,
    nx = nx,
    ny = 1)

parameters = {
    'time step duration': 10000,
    'substitutional molar volume': 1,
    'phase': {
	'name': "xi",
	'mobility': 1.,
	'gradient energy': 0.1,
	'value': 1.
    },
    'solvent': {
	'standard potential': Numeric.log(.7/.3),
	'barrier height': 1.
    }
}

parameters['substitutionals'] = (
    {
	'name': "c1",
	'diffusivity': 1.,
	'standard potential': Numeric.log(.3/.7),
	'barrier height': parameters['solvent']['barrier height'], 
	'value': 0.5
    },
)

fields = elphf.makeFields(mesh = mesh, parameters = parameters)

setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2)
fields['phase'].setValue(1.)
fields['phase'].setValue(0.,setCells)

equations = elphf.makeEquations(
    mesh = mesh, 
    fields = fields, 
    parameters = parameters
)

it = Iterator(equations = equations)

if __name__ == '__main__':
    phaseViewer = Gist1DViewer(vars = (fields['phase'],))
    concViewer = Gist1DViewer(vars = (fields['substitutionals'][0],))

    phaseViewer.plot()
    concViewer.plot()
	
    raw_input("press <return> to start...")

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    for i in range(50):
	it.timestep(1)
	
	phaseViewer.plot()
	concViewer.plot()

    # profile.stop()
	    
    raw_input("finished")

