#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 12/22/03 {4:28:58 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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
from equations.stdyConvDiffScEquation import SteadyConvectionDiffusionScEquation
from solvers.linearCGSSolver import LinearCGSSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from variables.cellVariable import CellVariable
from terms.powerLawConvectionTerm import PowerLawConvectionTerm
from terms.scSourceTerm import ScSourceTerm
from viewers.grid2DGistViewer import Grid2DGistViewer

d
    valueLeft = 0.
    valueRight = 1.

    L = parameters['L']
    nx = parameters['nx']
    ny = parameters['ny']

    mesh = Grid2D(L / nx, L / ny, nx, ny)

    var = CellVariable(
        name = "concentration",
        mesh = mesh,
        value = valueLeft)

    eq = SteadyConvectionDiffusionScEquation(
        var = var,
        diffusionCoeff = parameters['diffusion coeff'],
        convectionCoeff = parameters['convection coeff'],
        sourceCoeff = parameters['source coeff'],
        solver = LinearCGSSolver(tolerance = 1.e-15, steps = 2000),
        convectionScheme = parameters['convection scheme'],
        boundaryConditions =
        (
        FixedValue(faces = mesh.getFacesLeft(),value = valueLeft),
        FixedValue(faces = mesh.getFacesRight(),value = valueRight),
        FixedFlux(faces = mesh.getFacesTop(),value = 0.),
        FixedFlux(faces = mesh.getFacesBottom(),value = 0.)
	)
	)

    it = Iterator((eq,))
    
    return {
        'steps' : 1,
        'var' : var,
        'it' : it,
        'mesh' : mesh,
        'convection coeff' : convectionCoeff,
        'diffusion coeff' : diffusionCoeff,
        'source coeff' : sourceCoeff
	}
	
def runProcedure(parameters):
    newParameters = getParameters(parameters)

    it = newParameters['it']
    var = newParameters['var']
    steps = newParameters['steps']
    
    viewer = Grid2DGistViewer(var)
    it.timestep(steps)
    viewer.plot()
    raw_input()
