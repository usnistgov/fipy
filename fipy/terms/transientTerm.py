"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "transientTerm.py"
                                   created: 11/12/03 {11:36:25 AM} 
                               last update: 12/8/03 {3:57:08 PM} 
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

from cellTerm import CellTerm
import Numeric

class TransientTerm(CellTerm):
    def __init__(self,tranCoeff,mesh):
	weight = {'b vector': 0, 'new value': 1, 'old value': 1}
	CellTerm.__init__(self,weight,mesh) 
	self.tranCoeff = tranCoeff
	    
    def calculateCoeffGeom(self,dt):
	self.coeff = self.tranCoeff * self.mesh.getCellVolumes() / dt
	

