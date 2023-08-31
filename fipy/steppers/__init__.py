"""(Obsolete) utilities for iterating time steps

Use `steppyngstounes <https://pages.nist.gov/steppyngstounes/en/latest>`_
instead.
"""
from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.steppers.stepper import Stepper
from fipy.steppers.pseudoRKQSStepper import PseudoRKQSStepper
from fipy.steppers.pidStepper import PIDStepper
from fipy.tools.numerix import L1norm

__all__ = ["L1error", "L2error", "LINFerror", "sweepMonotonic",
           "error", "residual"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

def L1error(var, matrix, RHSvector):
    r"""
    .. math::

       \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_1}
       {\|\mathtt{var}^\text{old}\|_1}

    where :math:`\|\vec{x}\|_1` is the :math:`L^1` norm of :math:`\vec{x}`.

    Parameters
    ----------
    var : ~fipy.variables.cellVariable.CellVariable
        The `CellVariable` in question.
    matrix
        *(ignored)*
    RHSvector
        *(ignored)*
    """
    from fipy.tools.numerix import L1norm
    return error(var, matrix, RHSvector, L1norm)

def L2error(var, matrix, RHSvector):
    r"""
    .. math::

       \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_2}
       {\|\mathtt{var}^\text{old}\|_2}

    where :math:`\|\vec{x}\|_2` is the :math:`L^2` norm of :math:`\vec{x}`.

    Parameters
    ----------
    var : ~fipy.variables.cellVariable.CellVariable
        The `CellVariable` in question.
    matrix
        *(ignored)*
    RHSvector
        *(ignored)*
    """
    from fipy.tools.numerix import L2norm
    return error(var, matrix, RHSvector, L2norm)

def LINFerror(var, matrix, RHSvector):
    r"""
    .. math::

       \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_\infty}
       {\|\mathtt{var}^\text{old}\|_\infty}

    where :math:`\|\vec{x}\|_\infty` is the :math:`L^\infty` norm of :math:`\vec{x}`.

    Parameters
    ----------
    var : ~fipy.variables.cellVariable.CellVariable
        The `CellVariable` in question.
    matrix
        *(ignored)*
    RHSvector
        *(ignored)*
    """
    from fipy.tools.numerix import LINFnorm
    return error(var, matrix, RHSvector, LINFnorm)

def sweepMonotonic(fn, *args, **kwargs):
    """
    Repeatedly calls :func:`fn(*args, **kwargs)` until the residual returned by
    :func:`fn` is no longer decreasing.

    Parameters
    ----------
    fn : function
        The function to call
    *args
    **kwargs

    Returns
    -------
    float
    """
    oldres = 1e20
    # res is 0 the first time
    fn(*args, **kwargs)
    res = fn(*args, **kwargs)
    while res < oldres:
        oldres = res
        res = fn(*args, **kwargs)

    return res

def residual(var, matrix, RHSvector):
    r"""
    Determines the residual for the current solution matrix and variable.

    .. math::

        \|\mathsf{L}\vec{x} - \vec{b}\|_\infty

    where :math:`\|\vec{\xi}\|_\infty` is the :math:`L^\infty` norm of :math:`\vec{\xi}`.

    Parameters
    ----------
    var : ~fipy.variables.cellVariable.CellVariable
        The `CellVariable` in question, *prior* to solution.
    matrix : ~fipy.matrices.sparseMatrix._SparseMatrix
        The coefficient matrix at this step/sweep
    RHSvector : ndarray
        The right hand side vector
    """
    from fipy.tools.numerix import array, LINFnorm

    Lx = matrix * array(var)
    return LINFnorm(Lx - RHSvector)

def error(var, matrix, RHSvector, norm=L1norm):
    r"""
    .. math::

       \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_?}
       {\|\mathtt{var}^\text{old}\|_?}

    where :math:`\|\vec{x}\|_?` is the normalization of :math:`\vec{x}` provided
    by `norm`.

    Parameters
    ----------
    var : ~fipy.variables.cellVariable.CellVariable
        The `CellVariable` in question.
    matrix
        *(ignored)*
    RHSvector
        *(ignored)*
    norm : function
        A function that will normalize its `array` argument and return
        a single number (default: :func:`~fipy.tools.numerix.L1norm`).
    """
    denom = norm(var.old)
    return norm(var - var.old) / (denom + (denom == 0))
