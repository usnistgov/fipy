#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "transientTerm.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.terms.cellTerm import CellTerm
from fipy.terms import AlternativeMethodInBaseClass

class TransientTerm(CellTerm):
    r"""
    The `TransientTerm` represents
    
    .. math::

       \int_V \frac{\partial (\rho \phi)}{\partial t} dV \simeq
       \frac{(\rho_{P} \phi_{P} - \rho_{P}^\text{old} \phi_P^\text{old}) V_P}{\Delta t}
       
    where :math:`\rho` is the `coeff` value.

    The following test case verifies that variable coefficients and
    old coefficient values work correctly. We will solve the
    following equation

    .. math::

       \frac{ \partial \phi^2 } { \partial t } = k.
       
    The analytic solution is given by
    
    .. math::
        
       \phi = \sqrt{ \phi_0^2 + k t },
       
    where :math:`\phi_0` is the initial value.

    >>> phi0 = 1.
    >>> k = 1.
    >>> dt = 1.
    >>> relaxationFactor = 1.5
    >>> steps = 2
    >>> sweeps = 8
    
    >>> from fipy.meshes import Grid1D
    >>> mesh = Grid1D(nx = 1)
    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(mesh = mesh, value = phi0, hasOld = 1)
    >>> from fipy.terms.transientTerm import TransientTerm
    >>> from fipy.terms.implicitSourceTerm import ImplicitSourceTerm

    Relaxation, given by `relaxationFactor`, is required for a
    converged solution.
       
    >>> eq = TransientTerm(var) == ImplicitSourceTerm(-relaxationFactor) \
    ...                            + var * relaxationFactor + k 

    A number of sweeps at each time step are required to let the
    relaxation take effect.
    
    >>> for step in range(steps):
    ...     var.updateOld()
    ...     for sweep in range(sweeps):
    ...         eq.solve(var, dt = dt)

    Compare the final result with the analytical solution.
    
    >>> from fipy.tools import numerix
    >>> print var.allclose(numerix.sqrt(k * dt * steps + phi0**2))
    1
    """

    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
	return {
	    'b vector':  0, 
	    'new value': 1, 
	    'old value': 1,
            'diagonal': 0
	}
	
    def _calcGeomCoeff(self, mesh):
	return self.coeff * mesh.cellVolumes

    def _getTransientGeomCoeff(self, var):
        """
        Test to ensure that _getTransientGeomCoeff is not returning None when a
        TransientTerm is defined.

    def getTransientCoeff(self, mesh):
        return self._getGeomCoeff(mesh)

    def getTransientCoeff(self, mesh):
        return self._getGeomCoeff(mesh)
        
        >>> from fipy import *
        >>> from fipy.matrices.pysparseMatrix import _PysparseMatrix
        >>> m = Grid1D(nx=1)
        >>> var = CellVariable(mesh=m)
        >>> eq = TransientTerm(1) == ImplicitSourceTerm(1)
        >>> print CellVariable(mesh=m, value=eq._getTransientGeomCoeff(var))
        [ 1.]
        >>> eq.cacheMatrix()
        >>> eq.solve(var)
        >>> print not isinstance(eq.matrix, _PysparseMatrix) \
                or eq.matrix.asTrilinosMeshMatrix().numpyArray[0,0] == 1.
        True
        
        >>> eq = TransientTerm(-1) == ImplicitSourceTerm(1)
        >>> print CellVariable(mesh=m, value=eq._getTransientGeomCoeff(var))
        [-1.]
        >>> eq.cacheMatrix()
        >>> eq.solve(var)
        >>> print not isinstance(eq.matrix, _PysparseMatrix) \
                or eq.matrix.asTrilinosMeshMatrix().numpyArray[0,0] == -2.
        True

        """
        if CellTerm._getTransientGeomCoeff(self, var) is not None:
            AlternativeMethodInBaseClass('_getTransientGeomCoeff()')
        if var is self.var or self.var is None:
            return self._getGeomCoeff(var.mesh)
        else:
            return None

    @property
    def _transientVars(self):
        if len(super(TransientTerm, self)._transientVars) != 0:
            AlternativeMethodInBaseClass('_getDiffusionGeomCoeff()')
        return self._vars
    
        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

