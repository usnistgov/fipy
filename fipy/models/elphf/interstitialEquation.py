#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "interstitialEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 2/20/04 {6:58:54 PM} 
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

from concentrationEquation import ConcentrationEquation

class InterstitialEquation(ConcentrationEquation):
    def getConvectionCoeff(self, Cj, fields, diffusivity = None):
	if diffusivity is None:
	    diffusivity = Cj.getDiffusivity()
	Cj.weightedDiffusivity = (diffusivity * (1. + Cj.getHarmonicFaceValue())).transpose()
	
	return ConcentrationEquation.getConvectionCoeff(self, Cj, fields, Cj.weightedDiffusivity)
