#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 1/16/04 {9:35:02 AM} 
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

from tools.dimensions.physicalField import PhysicalField

import elphf

nx = 100
dx = "0.03 nm"
# L = nx * dx

mesh = Grid2D(
    dx = dx,
    dy = "1. m",
    nx = nx,
    ny = 1)
    
parameters = {
    'time step duration': "1e-12 s",
    'substitutional molar volume': "1.8e-5 m**3/mol",
    'phase': {
	'name': "xi",
	'mobility': "1e-2 m**3/J/s",
	'gradient energy': "3.6e-11 J/m",
	'value': 1.
    },
    'potential': {
	'name': "psi",
	'dielectric': 78.49
    },
    'solvent': {
	'standard potential': "34139.7265625 J/mol",
	'barrier height': "3.6e5 J/mol",
	'valence': 0
    }
}

parameters['interstitials'] = (
    {
	'name': "e-",
	'diffusivity': "1e-9 m**2/s",
	'standard potential': "-33225.9453125 J/mol",
	'barrier height': "0. J/mol",
	'valence': -1,
	'value': "111.110723815414 mol/l",
    },
)

parameters['substitutionals'] = (
    {
	'name': "SO4",
	'diffusivity': "1e-9 m**2/s",
	'standard potential': "24276.6640625 J/mol",
	'barrier height': parameters['solvent']['barrier height'],
	'valence': -2,
	'value': "0.000010414586295976 mol/l",
    },
    {
	'name': "Cu",
	'diffusivity': "1e-9 m**2/s",
	'standard potential': "-7231.81396484375 J/mol",
	'barrier height': parameters['solvent']['barrier height'],
	'valence': +2,
	'value': "55.5553718417909 mol/l",
    }
)

fields = elphf.makeFields(mesh = mesh, parameters = parameters)

setCells = mesh.getCells(lambda cell: cell.getCenter()[0] > mesh.getPhysicalShape()[0]/4.)
fields['phase'].setValue(0.,setCells)
fields['interstitials'][0].setValue("0.000111111503177394 MOLARVOLUME*mol/l", setCells)
fields['substitutionals'][0].setValue("0.249999982581341 MOLARVOLUME*mol/l", setCells)
fields['substitutionals'][1].setValue("0.249944439430068 MOLARVOLUME*mol/l", setCells)

it = elphf.makeIterator(mesh = mesh, fields = fields, parameters = parameters)

viewers = [Grid2DGistViewer(var = field) for field in fields['all']]

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
    it.timestep(steps = 1, maxSweeps = 5)
    
    for viewer in viewers:
	viewer.plot()

    raw_input()
    
# profile.stop()
	
raw_input()

