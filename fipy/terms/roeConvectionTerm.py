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

__all__ = ["RoeConvectionTerm"]

def MClimiter(theta):
    """
    Monotonized center limiter.
    """
    return numerix.maximum(0, numerix.minimum((1 + theta) / 2, numerix.minimum(2, 2 * theta)))



class _FirstOrderRoeVariable(FaceVariable):
    def __init__(self, var, coeff):
        super(FaceVariable, self).__init__(mesh=var.mesh, elementshape=(2,) + var.shape[:-1], cached=True)
        self.var = var
        self.coeff = self._requires(coeff)

    def _calcValue(self):
        from fipy.tools import smallMatrixVectorOps as smv
        
        eigenvalues, R, coeffDown, coeffUp, Ashape, id1, id2 = self._calcEigenvalues_()

        ## first order roe
            
        E = abs(eigenvalues)[:,:,numerix.newaxis] * numerix.identity(eigenvalues.shape[1])

        Abar = smv.mulinv(smv.mul(R, E), R)

        ## value.shape = (2, Nfac, Nequ, Nequ), first order 
        value = numerix.zeros((2,) + Ashape, 'd')

        value[0] = (coeffDown + Abar) / 2
        value[1] = (coeffUp - Abar) / 2

        self.R = R
        self.eigenvalues = eigenvalues
        self.E = E

        return value.transpose(0, 2, 3, 1)

    def _calcEigenvalues_(self):
        ## value[0] is the down implicit value
        ## value[1] is the up implicit value
        ## Imported when used to avoid calling cython unnecessarily.

        from fipy.tools import smallMatrixVectorOps as smv
        
        id1, id2 = self.mesh._adjacentCellIDs
        
        if len(self.coeff.shape) == 2:
            coeff = numerix.array(self.coeff)[:,numerix.newaxis,numerix.newaxis]
        else:
            coeff = numerix.array(self.coeff)
            
        ## coeff.shape = (Ndim, Nequ, Nequ, Ncell)

        coeffDown = (numerix.take(coeff, id1, axis=-1) * self.mesh._orientedAreaProjections[:, numerix.newaxis, numerix.newaxis]).sum(0).transpose(2, 0, 1)
        coeffUp = (numerix.take(coeff, id2, axis=-1) * self.mesh._orientedAreaProjections[:, numerix.newaxis, numerix.newaxis]).sum(0).transpose(2, 0, 1)

        ## A.shape = (Nfac, Nequ, Nequ)
        A = (coeffUp + coeffDown) / 2.

        eigenvalues, R = smv.sortedeig(A)

        return eigenvalues, R, coeffDown, coeffUp, A.shape, id1, id2

    @property
    def maxeigenvalue(self):
        eigenvalues, R, coeffDown, coeffUp, Ashape, id1, id2 = self._calcEigenvalues_()
        return  max(abs(eigenvalues).flat)

class _RoeVariable(FaceVariable):
    def __init__(self, var, coeff):
        super(FaceVariable, self).__init__(mesh=var.mesh, elementshape=(2,) + var.shape[:-1], cached=True)
        self.var = self._requires(var)
        self.firstOrderRoeVariable = _FirstOrderRoeVariable(var, coeff)
        
    @property
    def maxeigenvalue(self):
        return self.firstOrderRoeVariable.maxeigenvalue
        
    def _calcValue(self):
        value = self.firstOrderRoeVariable.value.transpose(0, 3, 1, 2).copy()

        ## second order correction with limiter
        correctionImplicit = self.getImplicitCorrection()
        
        value[0] -= correctionImplicit
        value[1] += correctionImplicit

        return value.transpose(0, 2, 3, 1)

    def getImplicitCorrection(self):
        R = self.firstOrderRoeVariable.R
        eigenvalues = self.firstOrderRoeVariable.eigenvalues
        E = self.firstOrderRoeVariable.E
        id1, id2 = self.mesh._adjacentCellIDs

        from fipy.tools import smallMatrixVectorOps as smv
        if len(self.var.shape) == 1:
            varT = numerix.array(self.var)[numerix.newaxis].transpose(1,0)
            varGradT = numerix.array(self.var.grad)[:,numerix.newaxis].transpose(2, 0, 1)
        else:
            varT = numerix.array(self.var).transpose(1,0)
            varGradT = numerix.array(self.var.grad).transpose(2, 0, 1)

        faceNormalsT = numerix.array(self.mesh._orientedFaceNormals).transpose(1,0)
        cellDistances = numerix.array(self.mesh._cellDistances)[:,numerix.newaxis]
        faceAreas = numerix.array(self.mesh._faceAreas)[:,numerix.newaxis]

        varDown = numerix.take(varT, id1, axis=0)
        varUp = numerix.take(varT, id2, axis=0)
        varDownGrad = smv.dot(numerix.take(varGradT, id1, axis=0), faceNormalsT[...,numerix.newaxis])
        varUpGrad = smv.dot(numerix.take(varGradT, id2, axis=0), faceNormalsT[...,numerix.newaxis])

        varUpUp = varDown + 2 * cellDistances * varUpGrad
        varDownDown = varUp - 2 * cellDistances * varDownGrad

##        varDiff = numerix.concatenate(((varUp - varDown)[...,numerix.newaxis],
##                                       (varDown - varDownDown)[...,numerix.newaxis],
##                                       (varUpUp - varUp)[...,numerix.newaxis]), axis=-1)



        alpha = smv.invmatvec(R, varUp - varDown)
        alphaDown = smv.invmatvec(R, varDown - varDownDown)
        alphaUp = smv.invmatvec(R, varUpUp - varUp)
        
        ##alpha = smv.invmul(R, varDiff[...,0][...,numerix.newaxis])
        ##alphaDown = smv.invmul(R, varDiff[...,1][...,numerix.newaxis])
        ##alphaUp = smv.invmul(R, varDiff[...,2][...,numerix.newaxis])

##         alphas = smv.invmul(R, varDiff)
##         alpha = alphas[...,0]
##         alphaDown = alphas[...,1]
##         alphaUp = alphas[...,2]

        alphaOther = numerix.where(eigenvalues > 0, alphaDown, alphaUp)

        theta = alphaOther / (alpha + 1e-20 * (alpha == 0))

        A = smv.mul(0.5 * E * (1 - float(self._dt) / (faceAreas * cellDistances)[...,numerix.newaxis] * E), R)
        
        return smv.mulinv(MClimiter(theta)[:,numerix.newaxis] * A, R)

    def _setdt(self, dt):
        if hasattr(self, '_dt') and not numerix.allclose(dt, self._dt):
            self._markStale()
        self._dt = dt

    def _getdt(self):
        return self._dt

    dt = property(_getdt, _setdt)

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

    def __init__(self, coeff=1.0, var=None):
        """
        :Parameters:
          - `coeff` : The `Term`'s coefficient value.

        """

        self.stencil = None
        
        if isinstance(coeff, _MeshVariable) and coeff.rank < 1:
            raise VectorCoeffError

        super(RoeConvectionTerm, self).__init__(coeff=coeff, var=var)

    def _getCoeffMatrix_(self, var, weight, dt):
        
        coeff = self._getGeomCoeff(var, dontCacheMe=False)
        coeff.dt = dt
        
        if self.coeffMatrix is None:
            self.coeffMatrix = {'cell 1 diag' : coeff[0],
                                'cell 1 offdiag': coeff[1],
                                'cell 2 diag': -coeff[1],
                                'cell 2 offdiag': -coeff[0]}

        return self.coeffMatrix

    def _getWeight(self, var, transientGeomCoeff, diffusionGeomCoeff):
        return {'implicit' : {'cell 1 diag' : 1}}
    
    def _calcGeomCoeff(self, var):
        return _RoeVariable(var, self.coeff)

    def maxeigenvalue(self, var):
        coeff = self._getGeomCoeff(var, dontCacheMe=False)
        coeff.dt = 1.
        return coeff.maxeigenvalue

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=1., transientGeomCoeff=None, diffusionGeomCoeff=None):

        var, L, b = super(RoeConvectionTerm, self)._buildMatrix(var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

        mesh = var.mesh
        if (not hasattr(self, 'constraintL')) or (not hasattr(self, 'constraintB')):

            constraintMask = var.faceGrad.constraintMask | var.faceValue.constraintMask

            diag =  self._getCoeffMatrix_(var, None, dt)['cell 1 diag'] * constraintMask
            offdiag = self._getCoeffMatrix_(var, None, dt)['cell 2 diag'] * constraintMask

            cellValue = _CellFaceValue(var)
            ghostValue = (cellValue + 2 * (var.faceValue - cellValue)) * constraintMask

            from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
            self.constraintL = _AddOverFacesVariable(diag) * mesh.cellVolumes
            self.constraintB = _AddOverFacesVariable(offdiag * ghostValue) * mesh.cellVolumes
        
        ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))
        L.addAt(numerix.array(self.constraintL).ravel(), ids.ravel(), ids.swapaxes(0,1).ravel())
        b += numerix.reshape(self.constraintB.value, ids.shape).sum(1).ravel()

        ## explicit
        
        b -= L * var.ravel()
        L = SparseMatrix(mesh=mesh)

        return (var, L, b)

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        from fipy.solvers import DefaultAsymmetricSolver
        solver = solver or super(RoeConvectionTerm, self)._getDefaultSolver(var, solver, *args, **kwargs)
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        return solver or DefaultAsymmetricSolver(*args, **kwargs)


