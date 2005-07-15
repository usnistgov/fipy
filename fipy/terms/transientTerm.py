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

    The following test case tests variable coefficients. We wish to solve the
    follwoing equation

    .. raw:: latex

        $$ \frac{ \partial \phi^2 } { \partial t } = k. $$ The analytic solution is given by
        $$ \phi = \sqrt{ \phi_0^2 + k t }, $$ where $\phi_0$

    is the initial value.

       >>> phi0 = 1.
       >>> k = 1.
       >>> dt = 1.
       
       >>> from fipy.meshes.grid1D import Grid1D
       >>> mesh = Grid1D(nx = 1)
       >>> from fipy.variables.cellVariable import CellVariable
       >>> var = CellVariable(mesh = mesh, value = phi0, hasOld = 1)
       >>> from fipy.terms.transientTerm import TransientTerm
       >>> eq = TransientTerm(var) - k

    We will do just one time step. Relaxation, given by `alpha`
    is required for a converged solution.
    
       >>> alpha = 0.5
       >>> for sweep in range(4):
       ...     tmpVar = var.copy()
       ...     eq.solve(var, dt = dt)
       ...     var.setValue(var * alpha + tmpVar * (1 - alpha))
       >>> import fipy.tools.numerix as numerix
       >>> print var.allclose(numerix.sqrt(k * dt + phi0**2))
       1
       
       
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
	
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

