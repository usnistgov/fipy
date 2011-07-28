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

##from fipy.terms.baseConvectionTerm import _BaseConvectionTerm
from fipy.variables.meshVariable import _MeshVariable
from fipy.terms.faceTerm import FaceTerm
from fipy.tools import numerix
from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.variables.faceVariable import FaceVariable
import scipy
import scipy.sparse
import scipy.sparse.linalg

class _RoeVariable(FaceVariable):
    def __init__(self, var, coeff, dt):
        super(FaceVariable, self).__init__(mesh=var.mesh, elementshape=(2,) + var.shape[:-1], cached=True)
        self.var = self._requires(var)
        self.dt = self._requires(dt)
        self.coeff = self._requires(coeff)
        
    def _calcValue(self):
        id1, id2 = self.mesh._adjacentCellIDs
        mesh = self.var.mesh
        
        coeffDown = (numerix.take(self.coeff, id1, axis=-1) * mesh._orientedAreaProjections[:, numerix.newaxis, numerix.newaxis]).sum(0)
        coeffUp = (numerix.take(self.coeff, id2, axis=-1) * mesh._orientedAreaProjections[:, numerix.newaxis, numerix.newaxis]).sum(0)
        
        ## Helper functions
        def mul(A, B):
            """
            Matrix multiply N, MxM matrices, A.shape = B.shape = (M, M, N).
            """
            return numerix.sum(A.swapaxes(0,1)[:, :, numerix.newaxis] * B[:, numerix.newaxis], 0)

        def inv(A):
            """
            Inverts N, MxM matrices, A.shape = (M, M, N).
            """
            return numerix.array(map(numerix.linalg.inv, A.transpose(2, 0, 1))).transpose(1, 2, 0)

        def mulinv(A, B):
            """
            Calculates the product (A.B^-1) of N, MxM matrices (A) and the inverse of N, MxM matrices (B).
            A.shape = (M, M, N), B.shape = (M, M, N).
            """
            tmp = B.transpose(2, 0, 1)
            N, M, M = tmp.shape
            import scipy
            BT = scipy.sparse.bsr_matrix((tmp, numerix.arange(N), numerix.arange(N + 1)), shape=(M * N, M * N)).transpose()
            PT = numerix.zeros((N, M, M), 'd')
            AT = A.transpose(2, 1, 0)
            
            for iCol in range(M):
                import scipy.sparse.linalg
                PT[:,:,iCol] = scipy.sparse.linalg.spsolve(BT.asformat("csc"), AT[:,:,iCol].ravel()).reshape((N, M))
            return PT.transpose(2, 1, 0)

        def invmul(A, B):
            """
            Calculates the product (A^-1.B) of N, MxM matrices (B) and the inverse of N, MxM matrices (A).
            A.shape = (M, M, N), B.shape = (M, M, N).
            """
            tmp = A.transpose(2, 0, 1)
            N, M, M = tmp.shape
            import scipy
            Amat = scipy.sparse.bsr_matrix((tmp, numerix.arange(N), numerix.arange(N + 1)), shape=(M * N, M * N))
            P = numerix.zeros((N, M, M), 'd')
            Bvec = B.transpose(2, 0, 1)
            
            for iCol in range(M):
                import scipy.sparse.linalg
                P[:,:,iCol] = scipy.sparse.linalg.spsolve(Amat.asformat("csc"), Bvec[:,:,iCol].ravel()).reshape((N, M))
            return P.transpose(1, 2, 0)

        def invmatvec(A, v):
            """
            Calculates the product (A^-1.v) of N, MxM matrices (B) and the inverse of N, MxM matrices (A).
            A.shape = (M, M, N), v.shape = (M, N).
            """
            M, M, N = A.shape
            tmp = A.transpose(2, 0, 1)
            import scipy
            Amat = scipy.sparse.bsr_matrix((tmp, numerix.arange(N), numerix.arange(N + 1)), shape=(M * N, M * N))
            P = scipy.sparse.linalg.spsolve(Amat.asformat("csc"), v.transpose(1, 0).ravel()).reshape((N, M))
            return P.transpose(1, 0)
            
        def eigvecN(A):
            """
            Calculate the eigenvalues and eigenvectors of N, MxM matrices, A.shape = (M, M, N).
            """
            tmp = zip(*map(numerix.linalg.eig, A.transpose(2, 0, 1)))
            return numerix.array(tmp[0]).swapaxes(0,1), numerix.array(tmp[1]).transpose(1,2,0)

        def vecN(A, eigs):
            """
            Calculates the eigenvectors for N, MxM matrices (A) using the eigenvalues (eigs).
            A.shape = (M, M, N), eigs.shape = (M, N).
            """
            M, M, N = A.shape
            vectors = numerix.zeros((N, M, M), 'd')
            for iCol in range(M):
                Amod = A - eigs[iCol] * numerix.identity(M)[..., numerix.newaxis]
                tmp = Amod.transpose(2, 0, 1)
                Amat = scipy.sparse.bsr_matrix((tmp, numerix.arange(N), numerix.arange(N + 1)), shape=(M * N, M * N))
                X = numerix.zeros((N, M), 'd')
                X[:,iCol] = 1.
                vectors[:,:,iCol] = scipy.sparse.linalg.gmres(Amat.asformat("csc"), numerix.zeros(M * N), X.ravel())[0].reshape((N, M))
                vectors[:,:,iCol] = vectors[:,:,iCol] / numerix.sqrt(numerix.sum(vectors[:,:,iCol]**2, 1))[...,numerix.newaxis]
            return vectors.transpose(1, 2, 0)

        def eigvec2(A):
            """
            Calculates the eigenvalues and eigenvectors of N, 2x2 matrices, A.shape = (2, 2, N).
            """
            tr = A[0,0] + A[1,1]
            det = A[0,0] * A[1,1] - A[1,0] * A[0,1]
            DeltaRt = numerix.sqrt(tr**2 - 4 * det)
            eigs = numerix.array(((tr + DeltaRt) / 2., (tr - DeltaRt) / 2.))
            return eigs, vecN(A, eigs)

        def eigvec(A):
            M, M, N = A.shape
            if M == 2:
                return eigvec2(A)
            else:
                return eigvecN(A)

        def sortedeig(A):
            """
            Caclulates the sorted eigenvalues and eigenvectors of N, MxM matrices, A.shape = (M, M, N).
            """
            N = A.shape[-1]
            
            eigenvalues, R = eigvec(A)

            order = eigenvalues.argsort(0).swapaxes(0, 1)
            Nlist = [[i] for i in xrange(N)]
            return (eigenvalues[order, Nlist].swapaxes(0, 1),
                    R[:, order, Nlist].swapaxes(1, 2))

        def matvec(A, v):
            """
            Multiply N, MxM matrices by N, M length vectors, A.shape = (M, M, N), v.shape = (M, N).
            """
            return numerix.sum(A * v, 1)

        def dot(v0, v1):
            """
            Dot product of two N, M-length vectors, v0.shape = (M, N), v1.shape = (M, N).
            """
            return numerix.sum(v0 * v1, 0)

        def MClimiter(theta):
            return numerix.maximum(0, numerix.minimum((1 + theta) / 2, 2, 2 * theta))
                               
        ## A.shape = (Nequ, Nequ, Nfac)
        A = (coeffUp + coeffDown) / 2.
        
        eigenvalues, R = sortedeig(A)        
        self._maxeigenvalue = max(abs(eigenvalues).flat)
        E = abs(eigenvalues) * numerix.identity(eigenvalues.shape[0])[..., numerix.newaxis]
        
        Abar = mulinv(mul(R, E), R)

        ## value.shape = (2, Nequ, Nequ, Nfac), first order 
        value = numerix.zeros((2,) + A.shape, 'd')
        value[0] = (coeffDown + Abar) / 2
        value[1] = (coeffUp - Abar) / 2

        ## second order correction with limiter
        varDown = numerix.take(self.var, id1, axis=-1)
        varUp = numerix.take(self.var, id2, axis=-1)
        varDownGrad = dot(numerix.take(self.var.grad(), id1, axis=-1), mesh._orientedFaceNormals)
        varUpGrad = dot(numerix.take(self.var.grad(), id2, axis=-1), mesh._orientedFaceNormals)
        
        varUpUp = varDown + 2 * mesh._cellDistances * varUpGrad
        varDownDown = varUp - 2 * mesh._cellDistances * varDownGrad

        alpha = invmatvec(R, varUp - varDown)
        alphaDown = invmatvec(R, varDown - varDownDown)
        alphaUp = invmatvec(R, varUpUp - varUp)
        
        alphaOther = numerix.where(eigenvalues > 0, alphaDown, alphaUp)
        theta = alphaOther / (alpha + 1e-20 * (alpha == 0))

        A = mul(0.5 * E * (1 - float(self.dt) / mesh._cellDistances * E), R)
        correctionImplicit = invmul(R.transpose(1, 0, 2), MClimiter(theta)[:, numerix.newaxis] * A.transpose(1, 0, 2)).transpose(1, 0, 2)

        value[0] -= correctionImplicit
        value[1] += correctionImplicit

        return value

    @property
    def maxeigenvalue(self):
        if not hasattr(self, '_maxeigenvalue'):
            self.value

        return self._maxeigenvalue

class _CellFaceValue(_CellToFaceVariable):
    def _calcValue_(self, alpha, id1, id2):
        return numerix.take(self.var, id1, axis=-1)

    @property
    def constraints(self):
        return super(_CellToFaceVariable, self).constraints

class RoeConvectionTerm(FaceTerm):
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

    def __init__(self, coeff=1.0, var=None, dt=1.0):
        """
        :Parameters:
          - `coeff` : The `Term`'s coefficient value.

        """

        self.stencil = None
        
        if isinstance(coeff, _MeshVariable) and coeff.rank < 1:
            raise VectorCoeffError

        super(RoeConvectionTerm, self).__init__(coeff=coeff, var=var)

        self.dt = dt
        
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
        return _RoeVariable(var, self.coeff, self.dt)

    def maxeigenvalue(self, var):
        return self._getGeomCoeff(var, dontCacheMe=False).maxeigenvalue

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=1., transientGeomCoeff=None, diffusionGeomCoeff=None):

        var, L, b = super(RoeConvectionTerm, self)._buildMatrix(var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

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
        ## explicit
##        b -= L * var.ravel()
##        L.matrix[:,:] = 0

        return (var, L, b)

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        solver = solver or super(RoeConvectionTerm, self)._getDefaultSolver(var, solver, *args, **kwargs)
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        return solver or DefaultAsymmetricSolver(*args, **kwargs)


