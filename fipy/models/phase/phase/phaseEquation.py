"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "phaseEquation.py"
                                   created: 11/12/03 {10:39:23 AM} 
                               last update: 11/26/03 {10:25:42 AM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
 Author: Daniel Wheeler
 E-mail: daniel.wheeler@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-12 JEG 1.0 original
###################################################################
"""

from matrixEquation import MatrixEquation
from terms.transientTerm import TransientTerm
from terms.explicitdiffusionTerm import ExplicitDiffusionTerm
from terms.scSourceTerm import ScSourceTerm
from terms.spSourceTerm import SpSourceTerm


class PhaseEquation(MatrixEquation):
    """
    Diffusion equation is implicit.
    """    
    def __init__(self,
                 var,
                 theta = 0.,
                 name='default_name',
                 solver='default_solver',
                 boundaryConditions=(),
                 parameters = {}):
        
        mesh = var.getMesh()
        self.parameters = parameters
        transientCoeff = self.parameters['tau']
        
	terms = (
	    TransientTerm(transientCoeff,mesh),
	    ExplicitDiffusionTerm(0.,mesh,boundaryConditions),
            SpSourceTerm(0.,mesh),
            ScSourceTerm(0.,mesh)
            )

	MatrixEquation.__init__(
            self,
            name,
            var,
            terms,
            solver)

    def preSolve():
        self.updateDiffCoeff()
        self.updateScCoeff()
        self.updateSpCoeff()

    def updateDiffCoeff():
        
