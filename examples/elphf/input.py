#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 12/22/03 {3:49:31 PM} 
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
import elphf
from meshes.grid2D import Grid2D
from componentVariable import ComponentVariable
from phaseVariable import PhaseVariable
# from variables.cellVariable import CellVariable
from viewers.grid2DGistViewer import Grid2DGistViewer

from profiler.profiler import Profiler
from profiler.profiler import calibrate_profiler

import Numeric

# valueLeft="0.3 mol/l"
# valueRight="0.4 mol/l"
# valueOther="0.2 mol/l"
valueLeft=0.3
valueRight=0.4
valueOther=0.2

nx = 40
dx = 1.
L = nx * dx

mesh = Grid2D(
    dx = dx,
    dy = 1.,
    nx = nx,
    ny = 1)
    
phase = PhaseVariable(
    name = "phase",
    mesh = mesh,
    value = 1.,
    viewer = Grid2DGistViewer
    )
    
var1 = ComponentVariable(
    name = "c1",
    standardPotential = Numeric.log(.3/.4) - Numeric.log(.1/.2),
    barrierHeight = 0.1,
    mesh = mesh,
    value = .35,
    viewer = Grid2DGistViewer
    )

var2 = ComponentVariable(
    name = "c2",
    standardPotential = Numeric.log(.4/.3) - Numeric.log(.1/.2),
    barrierHeight = 0.1,
    mesh = mesh,
    value = .35,
    viewer = Grid2DGistViewer
    )
   
var3 = ComponentVariable(
    name = "c3",
    standardPotential = Numeric.log(.2/.1) - Numeric.log(.1/.2),
    barrierHeight = 0.1,
    mesh = mesh,
    value = .15,
    viewer = Grid2DGistViewer
    )
   
rightCells = mesh.getCells(lambda center: center[0] > L/2.)

phase.setValue(0.,rightCells)
# var1.setValue(valueRight,rightCells)
# var2.setValue(valueLeft,rightCells)
# var3.setValue(0.1,rightCells)
    
fields = {
    'phase': phase,
    'substitutionals': (var1,var2,var3),
}

parameters = {
    'diffusivity': 1.,
    'time step duration': 10000.
}

it = elphf.makeIterator(mesh = mesh, fields = fields, parameters = parameters)

# print var1[:]
# print var2[:]
# print var3[:]

var1.plot()
var2.plot()
var3.plot()

raw_input()

# fudge = calibrate_profiler(10000)
# profile = Profiler('profile', fudge=fudge)

for i in range(50):
    it.iterate(1)
#     it.iterate(1,10000.)
    
#     print var1.getValue()
#     print var2.getValue()
    
#     var1.plot()
#     var2.plot()
#     var3.plot()

print var1
print var2
print var3
# profile.stop()
	
raw_input()

