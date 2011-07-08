#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "roeConvectionTerm.py"
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

from fipy.terms.baseConvectionTerm import _BaseConvectionTerm
from fipy.variables.meshVariable import _MeshVariable
from fipy.terms.faceTerm import FaceTerm

class RoeConvectionTerm(_BaseConvectionTerm):
    r"""
    A convection term implementing the Roe approximate Riemann flux update of the form

    .. math::

        F_i = \frac{1}{2} \left[ n_j u^-_{jik} \phi_k^- + n_j u^+_{jik} \phi_k^+ \right] -
        \frac{1}{2} n_j |\hat{A}_{jik}| \left( \phi_k^+ - \phi_k^- \right)

    where the flux is across a face between cell $C^-$ and cell $C^+$ and
    the normal $n_j$ points from $C^-$ to $C^+$. The question when
    implementing Roe solvers is how to approximate $\hat{A}_{jik}$ from the
    local Riemann problem

    .. math::

        \partial_t \phi_i + A_{jik} \partial_j \phi_k = 0

    where

    .. math::

        A_{jik} = \frac{ \partial \left( u_{jil} \phi_l \right) }{ \partial \phi_k }

    """

    def __init__(self, coeff=1.0, var=None):
        """
        :Parameters:
          - `coeff` : The `Term`'s coefficient value.

        """

        self.stencil = None
        
        if isinstance(coeff, _MeshVariable) and coeff.rank < 1:
            raise VectorCoeffError

        FaceTerm.__init__(self, coeff=coeff, var=var)
        
    def _getCoeffMatrix_(self, var, weight):
        if self.coeffMatrix is None:
            coeff = self._getGeomCoeff(var, dontCacheMe=False)
            self.coeffMatrix = {'cell 1 diag' : coeff[0],
                                'cell 1 offdiag': coeff[1],
                                'cell 2 diag': -coeff[1],
                                'cell 2 offdiag': -coeff[0]}
        return self.coeffMatrix

    def _getWeight(self, var, transientGeomCoeff, diffusionGeomCoeff):
        return {'implicit' : {'cell 1 diag' : 1}}
    
    def _calcGeomCoeff(self, var):
        from fipy.variables.roeVariable import _RoeVariable        
        return _RoeVariable(var, self.coeff)

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=1., transientGeomCoeff=None, diffusionGeomCoeff=None):
        from fipy.terms.faceTerm import FaceTerm
        return FaceTerm._buildMatrix(self, var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        solver = solver or super(RoeConvectionTerm, self)._getDefaultSolver(var, solver, *args, **kwargs)
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        return solver or DefaultAsymmetricSolver(*args, **kwargs)


