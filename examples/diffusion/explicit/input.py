#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 1/16/04 {10:52:51 AM} 
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

from meshes.grid2D import Grid2D
from equations.explicitDiffusionEquation import ExplicitDiffusionEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from variables.cellVariable import CellVariable
from viewers.grid2DGistViewer import Grid2DGistViewer

def getParameters(nx ,ny):

    valueLeft = 0.
    valueRight = 1.
    timeStep = 0.2

    mesh = Grid2D(1.,1.,nx,ny)

    var = CellVariable(
        name = "concentration",
        mesh = mesh,
        value = valueLeft)

    eq = ExplicitDiffusionEquation(
        var,
        transientCoeff = 1. / timeStep, 
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
        'timeStep' : timeStep,
        'tolerance' : 1e-3,
        'steps' : 10000,
        'mesh' : mesh,
        'var' : var,
        'it' : it
        }
        
    return parameters

if __name__ == '__main__':
    parameters = getParameters(10,1)
    it = parameters['it']
    var = parameters['var']
    steps = parameters['steps']
    viewer = Grid2DGistViewer(var)
    it.timestep(steps)
    viewer.plot()
    raw_input()

