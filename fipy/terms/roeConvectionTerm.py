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
from fipy.tools import numerix
from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.variables.faceVariable import FaceVariable

class _RoeVariable(FaceVariable):
    def __init__(self, var, coeff):
        super(FaceVariable, self).__init__(mesh=var.mesh, elementshape=(2,) + var.shape[:-1], cached=True)
        self.var = var

        self.coeff = self._requires(coeff)
        
    def _calcValue(self):
        id1, id2 = self.mesh._adjacentCellIDs
        mesh = self.var.mesh
        
        coeffDown = (numerix.take(self.coeff, id1, axis=-1) * mesh._orientedAreaProjections[:, numerix.newaxis, numerix.newaxis]).sum(0)
        coeffUp = (numerix.take(self.coeff, id2, axis=-1) * mesh._orientedAreaProjections[:, numerix.newaxis, numerix.newaxis]).sum(0)
        
        ## Helper functions
        def mul(A, B):
            """
            Matrix multiply N MxM matrices, A.shape = B.shape = (M, M, N).
            """
            return numerix.sum(A.swapaxes(0,1)[:, :, numerix.newaxis] * B[:, numerix.newaxis], 0)

        def inv(A):
            """
            Inverts N MxM matrices, A.shape = (M, M, N).
            """
            return numerix.array(map(numerix.linalg.inv, A.transpose(2, 0, 1))).transpose(1, 2, 0)

        def eig(A):
            """
            Calculate the eigenvalues and eigenvectors of N MxM matrices, A.shape = (M, M, N).
            """
            tmp = zip(*map(numerix.linalg.eig, A.transpose(2, 0, 1)))
            return numerix.array(tmp[0]).swapaxes(0,1), numerix.array(tmp[1]).transpose(1,2,0)

        def sortedeig(A):
            """
            Caclulates the sorted eigenvalues and eigenvectors of N MxM matrices, A.shape = (M, M, N).
            """
            N = A.shape[-1]
            eigenvalues, R = eig(A)
            order = eigenvalues.argsort(0).swapaxes(0, 1)
            Nlist = [[i] for i in xrange(N)]
            return (eigenvalues[order, Nlist].swapaxes(0, 1),
                    R[:, order, Nlist].swapaxes(1, 2))

        ## A.shape = (Nequ, Nequ, Nfac)
        A = (coeffUp + coeffDown) / 2.
        
        eigenvalues, R = sortedeig(A)
        E = abs(eigenvalues) * numerix.identity(eigenvalues.shape[0])[..., numerix.newaxis]
        Abar = mul(mul(R, E), inv(R))

        ## value.shape = (2, Nequ, Nequ, Nfac)+
        value = numerix.zeros((2,) + A.shape, 'd')
        value[0] = (coeffDown + Abar) / 2
        value[1] = (coeffUp - Abar) / 2

        return value

class _CellFaceValue(_CellToFaceVariable):
    def _calcValue_(self, alpha, id1, id2):
        return numerix.take(self.var, id1, axis=-1)

    @property
    def constraints(self):
        return super(_CellToFaceVariable, self).constraints

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
        return _RoeVariable(var, self.coeff)

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=1., transientGeomCoeff=None, diffusionGeomCoeff=None):
        from fipy.terms.faceTerm import FaceTerm
        var, L, b = FaceTerm._buildMatrix(self, var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

        mesh = var.mesh

        if (not hasattr(self, 'constraintL')) or (not hasattr(self, 'constraintB')):

            constraintMask = var.faceGrad.constraintMask | var.faceValue.constraintMask

            diag =  self._getCoeffMatrix_(var, None)['cell 1 diag'] * constraintMask
            offdiag = self._getCoeffMatrix_(var, None)['cell 2 diag'] * constraintMask

            cellValue = _CellFaceValue(var)
            ghostValue = (cellValue + 2 * (var.faceValue - cellValue)) * constraintMask

            from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
            self.constraintL = _AddOverFacesVariable(diag) * mesh.cellVolumes
            self.constraintB = _AddOverFacesVariable(offdiag * ghostValue) * mesh.cellVolumes

        ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))
        L.addAt(numerix.array(self.constraintL).ravel(), ids.ravel(), ids.swapaxes(0,1).ravel())
        b += numerix.reshape(self.constraintB.value, ids.shape).sum(1).ravel()
##        b -= L * var.ravel()
##        L.matrix[:,:] = 0
        return (var, L, b)

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        solver = solver or super(RoeConvectionTerm, self)._getDefaultSolver(var, solver, *args, **kwargs)
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        return solver or DefaultAsymmetricSolver(*args, **kwargs)


