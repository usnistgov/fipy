#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 12/29/03 {2:41:05 PM} 
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

"""Electrochemical Phase Field input file

    Build a mesh, variable, and diffusion equation with fixed (zero) flux
    boundary conditions at the top and bottom and fixed value boundary
    conditions at the left and right.
    
    Iterates a solution and plots the result with gist.
    
    Iteration is profiled for performance.
"""

import Numeric

from profiler.profiler import Profiler
from profiler.profiler import calibrate_profiler

from meshes.grid2D import Grid2D
from viewers.grid2DGistViewer import Grid2DGistViewer

import elphf

# valueLeft="0.3 mol/l"
# valueRight="0.4 mol/l"
# valueOther="0.2 mol/l"

nx = 1000
dx = 0.01
L = nx * dx

mesh = Grid2D(
    dx = dx,
    dy = 1.,
    nx = nx,
    ny = 1)
    
rightFunc = lambda cell: cell.getCenter()[0] > L/2.
    
parameters = {
    'diffusivity': 1.,
    'time step duration': 10000000.,
    'phase': {
	'name': "xi",
	'mobility': 1.,
	'gradient energy': 0.025,
	'initial': (
	    1.,
	    {
		'value': 0.,
		'func': rightFunc
	    }
	)
    },
    'potential': {
	'name': "psi",
	'permittivity': 1.5e1
    },
    'solvent': {
	'standard potential': Numeric.log(.4/.6),
	'barrier height': 10.0,
	'valence': +1
    }
}

parameters['interstitials'] = (
    {
	'name': "c1",
	'standard potential': Numeric.log(.1/.9),
	'barrier height': 0., #parameters['solvent']['barrier height'],
	'valence': -1,
	'initial': (0.5,),
    },
)

parameters['substitutionals'] = (
#     {
# 	'name': "c1",
# 	'standard potential': Numeric.log(.3/.4),
# 	'barrier height': parameters['solvent']['barrier height'],
# 	'initial': (0.35,),
#     },
    {
	'name': "c2",
	'standard potential': Numeric.log(.4/.3),
	'barrier height': parameters['solvent']['barrier height'],
	'valence': -1,
	'initial': (0.35,),
    },
    {
	'name': "c3",
	'standard potential': Numeric.log(.2/.1),
	'barrier height': parameters['solvent']['barrier height'],
	'valence': +1,
	'initial': (0.15,),
    }
)

fields = elphf.makeFields(mesh = mesh, parameters = parameters)

it = elphf.makeIterator(mesh = mesh, fields = fields, parameters = parameters)

viewers = [Grid2DGistViewer(field) for field in fields['all']]

for viewer in viewers:
    viewer.plot()
    
raw_input()

# fudge = calibrate_profiler(10000)
# profile = Profiler('profile', fudge=fudge)

# it.timestep(50)
# 
# for timestep in range(5):
#     it.sweep(50)
#     it.advanceTimeStep


for i in range(50):
    it.timestep(1)
#     raw_input()
    
    for viewer in viewers:
	viewer.plot()
    
# profile.stop()
	
raw_input()

