#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "stdyConvDiffEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 1/16/04 {11:28:44 AM} 
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
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

from fivol.equations.matrixEquation import MatrixEquation
from fivol.terms.transientTerm import TransientTerm
from fivol.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fivol.terms.powerLawConvectionTerm import PowerLawConvectionTerm

class SteadyConvectionDiffusionEquation(MatrixEquation):
    """
    Diffusion equation is implicit.
    """    
    def __init__(self,
                 var,
                 diffusionCoeff = 1.,
		 convectionCoeff = 1.,
                 solver='default_solver',
		 convectionScheme = PowerLawConvectionTerm,
                 boundaryConditions=()):
		     
	mesh = var.getMesh()
	
	diffusionTerm = ImplicitDiffusionTerm(diffusionCoeff,mesh,boundaryConditions)
	convectionTerm = convectionScheme(convectionCoeff, mesh, boundaryConditions, diffusionTerm)
		     
	terms = (
	    diffusionTerm,
	    convectionTerm
            )
	    
	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)

