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

from equations.matrixEquation import MatrixEquation
from terms.transientTerm import TransientTerm
from terms.explicitDiffusionTerm import ExplicitDiffusionTerm
from terms.scSourceTerm import ScSourceTerm
from terms.spSourceTerm import SpSourceTerm
from variables.variable import Variable
import Numeric

class PhaseEquation(MatrixEquation):
    """
    Diffusion equation is implicit.
    """    
    def __init__(self,
                 var,
                 theta = 0.,
                 temperature =0.,
                 solver = 'default_solver',
                 boundaryConditions = (),
                 parameters = {}):
        
        mesh = var.getMesh()
        if type(temperature) != type(var):
            temperature = Variable('temperature', mesh, temperature, hasOld = 0)
        self.temperature = temperature
        if type(theta) != type(var):
            theta = Variable('theta', mesh, theta, hasOld = 0)
        self.theta = theta
        self.parameters = parameters
        transientCoeff = self.parameters['tau']
        self.diffTerm = ExplicitDiffusionTerm(0.,mesh,boundaryConditions)
        self.spTerm = SpSourceTerm(0.,mesh)
        self.scTerm = ScSourceTerm(0.,mesh)

	terms = (
	    TransientTerm(transientCoeff,mesh),
	    self.diffTerm,
            self.scTerm,
            self.spTerm
            )

	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)

    def preSolve(self):
        self.updateDiffusion()
        self.updateSource()

    def updateDiffusion(self):
        alpha = self.parameters['alpha']
        N = self.parameters['symmetry']
        c2 = self.parameters['anisotropy']

        dphi = self.var.getFaceGradient()
        self.dphi = dphi
        z = Numeric.arctan2(dphi[:,1],dphi[:,0])
        z = N * (z - self.theta.getFaceArray())
        z = Numeric.tan(z / 2.)
        z = z * z;
        z = (1. - z) / (1. + z);
        z = (1.+ c2 * z);
        diffCoeff = alpha**2 * z * z;
        self.diffTerm.setDiffCoeff(diffCoeff)

    def updateSource(self):
        phi = self.var.getArray()
        t = self.temperature.getArray()

        ## driving force double well

        tmp = phi * (1 - phi)
        m = phi - 0.5 + t * tmp

        sc = (m > 0.) * m * phi
        sp = m * (phi - (m < 0.))
    
        ## theta source terms

        thetaMag = self.theta.getOld().getGradientMagnitude()
        s = self.parameters['s']
        epsilon = self.parameters['epsilon']

        sp += (2*s + epsilon**2 * thetaMag) * thetaMag

        ## anisotropy

##        z = Numeric.atan2(self.dphi[:,1],self.dphi[:,0]);
##        z = N * (z-self.theta);
##        z = tan(0.5 * z);
##        zsq = z * z;
##        b = (1-zsq) / (1+zsq);
##        db = -N * 2. *z / (1+zsq);
##        ff = alphasq * c2 * (1.+c2 * b) * db
        
##        sc + = phaseTools.add_over_faces_inline(self.ff,-self.dphi[:,1],self.dphi[:,0],mesh)

        self.scTerm.setScCoeff(sc)
        self.spTerm.setSpCoeff(sp)

