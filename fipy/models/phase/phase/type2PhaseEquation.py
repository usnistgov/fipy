"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "type1PhaseEquation.py"
                                   created: 12/24/03 {10:39:23 AM} 
                               last update: 12/24/03 {5:53:13 PM} 
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

from phaseEquation import PhaseEquation
import Numeric

class Type2PhaseEquation(PhaseEquation):
    """
    Diffusion equation is implicit.
    """    
    def __init__(self,
                 var,
                 solver = 'default_solver',
                 boundaryConditions = (),
                 fields = {},
                 parameters = {}):

	temp = fields['temperature']
        kappa1 = parameters['kappa 1']
        kappa2 = parameters['kappa 2']
        pi = Numeric.pi
        
        mPhi = var - 0.5 - kappa1 / pi * Numeric.atan(kappa2 * temp)

        PhaseEquation.__init__(self,
                               var = var,
                               solver = solver,
                               boundaryConditions = boundaryConditions,
                               fields = fields,
                               parameters = parameters,
                               mPhi = mPhi)
        
