#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 7/26/04 {8:34:24 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
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

import Numeric

from fipy.tools.profiler.profiler import Profiler
from fipy.tools.profiler.profiler import calibrate_profiler

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
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
    'substitutional molar volume': 1.,
    'phase': {
	'name': "xi",
	'mobility': 1.,
	'gradient energy': 0.025,
	'value': 1.
    },
    'solvent': {
	'standard potential': Numeric.log(.4/.6) - Numeric.log(1.3/1.4),
	'barrier height': 1.
    }
}

parameters['interstitials'] = (
    {
	'name': "c1",
	'diffusivity': 1.,
	'standard potential': Numeric.log(.3/.4) - Numeric.log(1.3/1.4),
	'barrier height': 0., 
	'value': 0.35
    },
)

parameters['substitutionals'] = (
    {
	'name': "c2",
	'diffusivity': 1.,
	'standard potential': Numeric.log(.4/.3) - Numeric.log(1.3/1.4),
	'barrier height': parameters['solvent']['barrier height'], 
	'value': 0.35
    },
    {
	'name': "c3",
	'diffusivity': 1.,
	'standard potential': Numeric.log(.2/.1) - Numeric.log(1.3/1.4),
	'barrier height': parameters['solvent']['barrier height'], 
	'value': 0.15
    },
)

fields = elphf.makeFields(mesh = mesh, parameters = parameters)

setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2.)
fields['phase'].setValue(1.)
fields['phase'].setValue(0.,setCells)

equations = elphf.makeEquations(
    mesh = mesh, 
    fields = fields, 
    parameters = parameters
)

it = Iterator(equations = equations)

if __name__ == '__main__':
    viewers = [Grid2DGistViewer(var = field) for field in fields['all']]

    for viewer in viewers:
	viewer.plot()
	
    raw_input()

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    for i in range(500):
	it.timestep(1)
	
	for viewer in viewers:
	    viewer.plot()
	
    # profile.stop()
	    
    raw_input()

