#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 12/16/03 {9:42:04 AM} 
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
from variables.cellVariable import CellVariable
from viewers.grid2DGistViewer import Grid2DGistViewer

from profiler.profiler import Profiler
from profiler.profiler import calibrate_profiler

valueLeft="0.3 mol/l"
valueRight="0.6 mol/l"
# valueLeft=0.3
# valueRight=0.6

nx = 40
dx = 1.
L = nx * dx

mesh = Grid2D(
    dx = dx,
    dy = 1.,
    nx = nx,
    ny = 40)

var1 = CellVariable(
    name = "c1",
    mesh = mesh,
    value = valueLeft,
    viewer = Grid2DGistViewer
    )

var2 = CellVariable(
    name = "c2",
    mesh = mesh,
    value = valueRight,
    viewer = Grid2DGistViewer
    )
   
rightCells = mesh.getCells(lambda center: center[0] > L/2.)

var1.setValue(valueRight,rightCells)
var2.setValue(valueLeft,rightCells)
    
parameters = {
    'substitutionals': (var1,var2),
    'diffusivity': 1.
}

it = elphf.makeIterator(mesh = mesh, parameters = parameters)

var1.plot()
var2.plot()

fudge = calibrate_profiler(10000)
profile = Profiler('profile', fudge=fudge)

for i in range(5):
    it.iterate(1,10000.)
    
#     print var1.getValue()
#     print var2.getValue()
    
    var1.plot()
    var2.plot()

profile.stop()
	
raw_input()

