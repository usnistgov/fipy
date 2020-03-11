from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.terms.abstractConvectionTerm import _AbstractConvectionTerm
from fipy.variables.faceVariable import FaceVariable
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import inline
from fipy.tools import numerix

class _UpwindConvectionTermAlpha(FaceVariable):
    def __init__(self, P):
        FaceVariable.__init__(self, mesh=P.mesh, elementshape=P.shape[:-1])
        self.P = self._requires(P)

    if inline.doInline:
        def _calcValue(self):
            P  = self.P.numericValue
            alpha = self._array.copy()

            inline._runInline("""
                alpha[i] = 0.5;

                if (P[i] > 0.) {
                    alpha[i] = 1.;
                } else {
                    alpha[i] = 0.;
                }
            """,
            alpha=alpha, P=P,
            ni = len(P.flat))

            return self._makeValue(value=alpha)
    else:
        def _calcValue(self):
            P  = self.P.numericValue
            alpha = numerix.where(P > 0., 1., 0.)
            return PhysicalField(value=alpha)

class _AbstractUpwindConvectionTerm(_AbstractConvectionTerm):
    def _alpha(self, P):
        return _UpwindConvectionTermAlpha(P)

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):

        var, L, b = _AbstractConvectionTerm._buildMatrix(self, var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

        mesh = var.mesh
        
        ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))
        
        if (hasattr(self, 'constraintL')) or (hasattr(self, 'constraintB')):
            # remove L from matrix and b from right hand side
            # necessary due to call of _buildMatrix on super class
            L.addAt(numerix.array(-self.constraintL).ravel(), ids.ravel(), ids.swapaxes(0, 1).ravel())
            b -= numerix.reshape(self.constraintB.value, ids.shape).sum(0).ravel()          
            
            constraintMask = var.faceGrad.constraintMask | var.arithmeticFaceValue.constraintMask
            
            weight = self._getWeight(var, transientGeomCoeff, diffusionGeomCoeff)

            if 'implicit' in weight:
                alpha = weight['implicit']['cell 1 diag']
            else:
                alpha = 0.0

            exteriorCoeff =  self.coeff * mesh.exteriorFaces
            
            
            self.constraintL = (alpha * constraintMask * exteriorCoeff).divergence * mesh.cellVolumes
            self.constraintB =  -((1 - alpha) * var.arithmeticFaceValue * constraintMask * exteriorCoeff).divergence * mesh.cellVolumes
        
        L.addAt(numerix.array(self.constraintL).ravel(), ids.ravel(), ids.swapaxes(0, 1).ravel())
        b += numerix.reshape(self.constraintB.value, ids.shape).sum(0).ravel()

        return (var, L, b)
      
def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
