#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input1DphaseQuaternary.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 7/29/04 {2:40:36 PM} 
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
 # protection and is in the public domain.  PFM is an experimental
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

"""
This example adds two more substitutional components to `input1DphaseBinary.py`.

We start with uniform concentration fields

.. raw:: latex

   $C_1 = C_2 = 0.35$ and $C_3 = 0.15$.
   
Again, this problem does not have an analytical solution, so after
iterating to equilibrium

    >>> for step in range(40):
    ...     it.timestep()

we confirm that the far-field phases have remained separated

    >>> Numeric.allclose(Numeric.take(fields['phase'], (0,-1)), (1.0, 0.0), rtol = 2e-3, atol = 2e-3)
    1
    
and that the concentration fields has appropriately segregated into into
their respective phases

    >>> Numeric.allclose(Numeric.take(fields['substitutionals'][0], (0,-1)), (0.4, 0.3), rtol = 2e-3, atol = 2e-3)
    1
    >>> Numeric.allclose(Numeric.take(fields['substitutionals'][1], (0,-1)), (0.3, 0.4), rtol = 2e-3, atol = 2e-3)
    1
    >>> Numeric.allclose(Numeric.take(fields['substitutionals'][2], (0,-1)), (0.1, 0.2), rtol = 2e-3, atol = 2e-3)
    1
"""
__docformat__ = 'restructuredtext'
 
import Numeric

from fipy.tools.profiler.profiler import Profiler
from fipy.tools.profiler.profiler import calibrate_profiler

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
	'gradient energy': 0.025,
	'value': 1.
    },
    'solvent': {
	'standard potential': Numeric.log(.1/.2),
	'barrier height': 1.
    }
}

parameters['substitutionals'] = (
    {
	'name': "c1",
	'diffusivity': 1.,
	'standard potential': Numeric.log(.3/.4),
	'barrier height': parameters['solvent']['barrier height'], 
	'value': 0.35
    },
    {
	'name': "c2",
	'diffusivity': 1.,
	'standard potential': Numeric.log(.4/.3),
	'barrier height': parameters['solvent']['barrier height'], 
	'value': 0.35
    },
    {
	'name': "c3",
	'diffusivity': 1.,
	'standard potential': Numeric.log(.2/.1),
	'barrier height': parameters['solvent']['barrier height'], 
	'value': 0.15
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
    concViewer = Gist1DViewer(vars = fields['substitutionals'])

    phaseViewer.plot()
    concViewer.plot()
	
    raw_input()

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    for i in range(50):
	it.timestep(1)
	
	phaseViewer.plot()
	concViewer.plot()
	
    # profile.stop()
	    
    raw_input()

