from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable

__all__ = ["ModularVariable"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ModularVariable(CellVariable):
    r"""
    The `ModularVariable` defines a variable that exists on the circle between
    :math:`-\pi` and :math:`\pi`

    The following examples show how `ModularVariable` works. When
    subtracting the answer wraps back around the circle.

    >>> from fipy.meshes import Grid1D
    >>> mesh = Grid1D(nx = 2)
    >>> from fipy.tools import numerix
    >>> pi = numerix.pi
    >>> v1 = ModularVariable(mesh = mesh, value = (2*pi/3, -2*pi/3))
    >>> v2 = ModularVariable(mesh = mesh, value = -2*pi/3)
    >>> print(numerix.allclose(v2 - v1, (2*pi/3, 0)))
    1

    Obtaining the arithmetic face value.

    >>> print(numerix.allclose(v1.arithmeticFaceValue, (2*pi/3, pi, -2*pi/3)))
    1

    Obtaining the gradient.

    >>> print(numerix.allclose(v1.grad, ((pi/3, pi/3),)))
    1

    Obtaining the gradient at the faces.

    >>> print(numerix.allclose(v1.faceGrad, ((0, 2*pi/3, 0),)))
    1

    Obtaining the gradient at the faces but without modular
    arithmetic.

    >>> print(numerix.allclose(v1.faceGradNoMod, ((0, -4*pi/3, 0),)))
    1
    """

    _modIn = """
    # define pi 3.141592653589793
    # define mod(x) (fmod(x + 3. * pi, 2. * pi) - pi)
    """

    def _setValueInternal(self, value, unit=None, array=None):
        """
        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(nx = 4)
        >>> from fipy.variables.modularVariable import ModularVariable
        >>> var = ModularVariable(mesh = mesh, value = 1., hasOld = 1)
        >>> answer = CellVariable(mesh=mesh, value=1.)
        >>> print(var.allclose(answer))
        True
        >>> var.value = 1
        >>> print(var.allclose(answer))
        True
        """
        value = self._makeValue(value=value, unit=unit, array=array)
        from fipy.variables.modPhysicalField import _ModPhysicalField
        self._value = _ModPhysicalField(value=value, unit=unit, array=array)

    def updateOld(self):
        """
        Set the values of the previous solution sweep to the current values.
        Test case due to bug.

        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(nx = 1)
        >>> var = ModularVariable(mesh=mesh, value=1., hasOld=1)
        >>> var.updateOld()
        >>> var[:] = 2
        >>> answer = CellVariable(mesh=mesh, value=1.)
        >>> print(var.old.allclose(answer))
        True
        """
        self.value = (self.value.mod(self().inRadians()))
        if self._old is not None:
            self._old.value = (self._value.value.copy())

    @property
    def grad(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `CellVariable` (first-order
        gradient). Adjusted for a `ModularVariable`
        """
        if not hasattr(self, '_grad'):
            from fipy.variables.modCellGradVariable import _ModCellGradVariable
            self._grad = _ModCellGradVariable(self, self._modIn, self._value.mod)

        return self._grad

    @property
    def arithmeticFaceValue(self):
        r"""
        Returns a `FaceVariable` whose value corresponds to the arithmetic interpolation
        of the adjacent cells:

        .. math::

           \phi_f = (\phi_1 - \phi_2) \frac{d_{f2}}{d_{12}} + \phi_2

        Adjusted for a `ModularVariable`
        """
        if not hasattr(self, '_arithmeticFaceValue'):
            from fipy.variables.modCellToFaceVariable import _ModCellToFaceVariable
            self._arithmeticFaceValue = _ModCellToFaceVariable(self, self._modIn)

        return self._arithmeticFaceValue

    @property
    def faceGrad(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `FaceVariable` (second-order
        gradient). Adjusted for a `ModularVariable`
        """
        if not hasattr(self, '_faceGrad'):
            from fipy.variables.modFaceGradVariable import _ModFaceGradVariable
            self._faceGrad = _ModFaceGradVariable(self, self._modIn)

        return self._faceGrad

    @property
    def faceGradNoMod(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `FaceVariable` (second-order
        gradient). Not adjusted for a `ModularVariable`
        """

        if not hasattr(self, '_faceGradNoMod'):
            class NonModularTheta(CellVariable):
                def __init__(self, modVar):
                    CellVariable.__init__(self, mesh = modVar.mesh)
                    self.modVar = self._requires(modVar)

                def _calcValue(self):
                    return self.modVar.value

            self._faceGradNoMod = NonModularTheta(self).faceGrad

        return self._faceGradNoMod

    def __sub__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return -other + self
        else:
            return self._BinaryOperatorVariable(lambda a, b: a-b, other, canInline=False)

    def __rsub__(self, other):
        return self._BinaryOperatorVariable(lambda a, b: b-a, other, canInline=False)

    def _getArithmeticBaseClass(self, other=None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            return ModularVariable

        return CellVariable._getArithmeticBaseClass(self, other)


def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

