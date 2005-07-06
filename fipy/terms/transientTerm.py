#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "transientTerm.py"
 #                                    created: 11/12/03 {11:36:25 AM} 
 #                                last update: 7/6/05 {4:26:44 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
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

__docformat__ = 'restructuredtext'

from fipy.terms.cellTerm import CellTerm

class TransientTerm(CellTerm):
    r"""
    The `TransientTerm` is discretized in the following way
    
    .. raw:: latex

       $$ \int_V \frac{\partial (\rho \phi)}{\partial t} dV \simeq
       \frac{\rho_{P}(\phi_{P} - \phi_P^\text{old}) V_P}{\Delta t} $$
       where $\rho$ is the

    `coeff` value.

    Usage ::

        TransientTerm(coeff = <CellVariable|Float>)
        
    """

    def _getWeight(self, mesh):
	return {
	    'b vector':  0, 
	    'new value': 1, 
	    'old value': 1,
            'diagonal': 0
	}
	
    def _calcGeomCoeff(self, mesh):
	self.geomCoeff = self.coeff * mesh.getCellVolumes()
	

