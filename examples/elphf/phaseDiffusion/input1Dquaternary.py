#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input1DphaseQuaternary.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 1/16/04 {11:41:23 AM} 
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

from fivol.profiler.profiler import Profiler
from fivol.profiler.profiler import calibrate_profiler

from fivol.meshes.grid2D import Grid2D
from fivol.viewers.grid2DGistViewer import Grid2DGistViewer

import fivol.examples.elphf.elphf as elphf

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

setCells = mesh.getCells(lambda cell: cell.getCenter()[0] > L/2)
fields['phase'].setValue(1.)
fields['phase'].setValue(0.,setCells)

it = elphf.makeIterator(mesh = mesh, fields = fields, parameters = parameters)

if __name__ == '__main__':
    viewers = [Grid2DGistViewer(var = field) for field in fields['all']]

    for viewer in viewers:
	viewer.plot()
	
    raw_input()

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    for i in range(50):
	it.timestep(1)
	
	for viewer in viewers:
	    viewer.plot()
	
    # profile.stop()
	    
    raw_input()

