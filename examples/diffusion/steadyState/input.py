#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 1/16/04 {11:58:39 AM} 
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

from fivol.profiler.profiler import Profiler
from fivol.profiler.profiler import calibrate_profiler
 
from fivol.meshes.grid2D import Grid2D
from fivol.equations.diffusionEquation import DiffusionEquation
from fivol.solvers.linearPCGSolver import LinearPCGSolver
from fivol.boundaryConditions.fixedValue import FixedValue
from fivol.boundaryConditions.fixedFlux import FixedFlux
from fivol.iterators.iterator import Iterator
from fivol.variables.cellVariable import CellVariable
from fivol.viewers.grid2DGistViewer import Grid2DGistViewer

def getParameters(nx ,ny):

    valueLeft = 0.
    valueRight = 1.

    mesh = Grid2D(1.,1.,nx,ny)

    var = CellVariable(
        name = "concentration",
        mesh = mesh,
        value = valueLeft)

    eq = DiffusionEquation(
        var,
        transientCoeff = 0., 
        diffusionCoeff = 1.,
        solver = LinearPCGSolver(
        tolerance = 1.e-15, 
        steps = 1000
        ),
        boundaryConditions=(
        FixedValue(mesh.getFacesLeft(),valueLeft),
        FixedValue(mesh.getFacesRight(),valueRight),
        FixedFlux(mesh.getFacesTop(),0.),
        FixedFlux(mesh.getFacesBottom(),0.)
        )
        )

    it = Iterator((eq,))

    parameters = {
        'valueLeft' : valueLeft,
        'valueRight' : valueRight,
        'tolerance' : 1e-8,
        'steps' : 1,
        'mesh' : mesh,
        'var' : var,
        'it' : it
        }
        
    return parameters

if __name__ == '__main__':
    parameters = getParameters(50,50)
    it = parameters['it']
    var = parameters['var']
    viewer = Grid2DGistViewer(var)
    
    fudge = calibrate_profiler(10000)
    profile = Profiler('profile', fudge=fudge)
    
    it.timestep(1)
    
    profile.stop()
    
    viewer.plot()
    raw_input()
