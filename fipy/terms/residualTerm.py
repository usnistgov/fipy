from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm
from fipy.variables.cellVariable import CellVariable

__all__ = ["ResidualTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ResidualTerm(_ExplicitSourceTerm):
    r"""

    The `ResidualTerm` is a special form of explicit `SourceTerm` that adds the
    residual of one equation to another equation. Useful for Newton's method.
    """
    def __init__(self, equation, underRelaxation=1.):
        self.equation = equation
        self.underRelaxation = underRelaxation

        _ExplicitSourceTerm.__init__(self, var=None)

    def __repr__(self):
        return r"$\Delta$[" + repr(self.equation) + "]"

    def _getGeomCoeff(self, var):
        return self.coeff

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        vec = self.equation.justResidualVector(var=None,
                                               boundaryConditions=boundaryConditions,
                                               dt=dt)

        self.coeff = CellVariable(mesh=var.mesh, value=vec * self.underRelaxation)
        self.geomCoeff = None
        self.coeffVectors = None

        return _ExplicitSourceTerm._buildMatrix(self, var=var, SparseMatrix=SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)
