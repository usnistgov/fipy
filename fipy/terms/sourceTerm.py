from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.cellTerm import CellTerm
from fipy.terms import AbstractBaseClassError
from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix

__all__ = ["SourceTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SourceTerm(CellTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff=0., var=None):
        r"""
        Parameters
        ----------
        coeff : float or ~fipy.variables.cellVariable.CellVariable
            Coefficient of source (default: 0)
        var : ~fipy.variables.cellVariable.CellVariable
            Variable :math:`\phi` that
            :class:`~fipy.terms.sourceTerm.SourceTerm` is implicit in.
        """
        if self.__class__ is SourceTerm:
            raise AbstractBaseClassError
        super(SourceTerm, self).__init__(coeff=coeff, var=var)

    def _calcGeomCoeff(self, var):
        self._checkCoeff(var)

        if self.coeff.shape != () and self.coeff.shape[-1] != len(var.mesh.cellVolumes):
            return self.coeff[..., numerix.newaxis] * CellVariable(mesh=var.mesh, value=var.mesh.cellVolumes)
        else:
            return self.coeff * CellVariable(mesh=var.mesh, value=var.mesh.cellVolumes)

    def _checkDt(self, dt):
        return 1.
