#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 5/5/04 {6:40:39 PM} 
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

from fipy.tools.profiler.profiler import Profiler
from fipy.tools.profiler.profiler import calibrate_profiler

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.iterators.iterator import Iterator

import fipy.models.elphf.elphf as elphf


nx = 400
dx = 0.01
L = nx * dx

parameters = {
    'mesh': {
	'nx': nx,
	'ny': 1,
	'dx': dx,
	'dy': dx
    },
    'time step duration': 10000,
    'substitutional molar volume': 1,
    'phase': {
	'name': "xi",
	'mobility': 1.,
	'gradient energy': 0.025,
	'value': 1.
    },
    'solvent': {
	'standard potential': 0.,
	'barrier height': 1.
    }
}

mesh = Grid2D(
    dx = parameters['mesh']['dx'],
    dy = parameters['mesh']['dy'],
    nx = parameters['mesh']['nx'],
    ny = parameters['mesh']['ny'])
    
fields = elphf.makeFields(mesh = mesh, parameters = parameters)

setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2)
fields['phase'].setValue(1.)
fields['phase'].setValue(0.,setCells)

equations, timeStepDuration = elphf.makeEquations(
    mesh = mesh, 
    fields = fields, 
    parameters = parameters
)

it = Iterator(equations = equations, timeStepDuration = timeStepDuration)

if __name__ == '__main__':
    viewers = [Grid2DGistViewer(var = field) for field in fields['all']]

    for viewer in viewers:
	viewer.plot()
	
    raw_input()

    fudge = calibrate_profiler(10000)
    profile = Profiler('profile', fudge=fudge)

    for i in range(50):
	it.timestep(1)
	
## 	for viewer in viewers:
## 	    viewer.plot()
	
    profile.stop()
	    
    raw_input()

