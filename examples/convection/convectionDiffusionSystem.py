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
from iterators.iterator import Iterator
from variables.cellVariable import CellVariable
from terms.powerLawConvectionTerm import PowerLawConvectionTerm
from terms.scSourceTerm import ScSourceTerm
from viewers.grid2DGistViewer import Grid2DGistViewer

class ConvectionDiffusionSystem:

    def __init__(self):

        self.steps = 1
        self.valueLeft = 0.
        self.valueRight = 1.

        self.mesh = Grid2D(self.L / self.nx, self.L / self.ny, self.nx, self.ny)
        
        self.var = CellVariable(
            name = "concentration",
            mesh = self.mesh,
            value = self.valueLeft)
        
        eq = SteadyConvectionDiffusionScEquation(
            var = self.var,
            diffusionCoeff = self.diffCoeff,
            convectionCoeff = self.convCoeff,
            sourceCoeff = self.sourceCoeff,
            solver = LinearCGSSolver(tolerance = 1.e-15, steps = 2000),
            convectionScheme = self.convectionScheme,
            boundaryConditions = self.getBoundaryConditions()
            )
        
        self.it = Iterator((eq,))

        self.parameters = {
            'steps' : 1,
            'var' : self.var,
            'it' : self.it,
            'mesh' : self.mesh,
            'convection coeff' : self.convCoeff,
            'diffusion coeff' : self.diffCoeff,
            'source coeff' : self.sourceCoeff,
            'L' : self.L
            }


    def getParameters(self):
        return self.parameters
    
    def run(self):
        viewer = Grid2DGistViewer(self.var)
        self.it.timestep(self.steps)
        viewer.plot()
        raw_input()
