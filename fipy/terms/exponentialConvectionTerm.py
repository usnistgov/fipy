#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "exponentialConvectionTerm.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.terms.asymmetricConvectionTerm import _AsymmetricConvectionTerm
from fipy.variables.faceVariable import FaceVariable

__all__ = ["ExponentialConvectionTerm"]

class _ExponentialConvectionTermAlpha(FaceVariable):
    def __init__(self, P):
        FaceVariable.__init__(self, P.mesh)
        self.P = self._requires(P)

    def _calcValue(self):
        """

        Test case added because `and` was being used instead of bitwise `&`.

            >>> from fipy.meshes import Grid1D
            >>> mesh = Grid1D(nx = 3)
            >>> from fipy.variables.faceVariable import FaceVariable
            >>> P = FaceVariable(mesh = mesh, value = (1e-3, 1e+71, 1e-3, 1e-3))
            >>> alpha = ExponentialConvectionTerm([1])._alpha(P)
            >>> print alpha
            [ 0.5  1.   0.5  0.5]

        """
        eps = 1e-3
        largeValue = 101.0
        P  = self.P.numericValue

        P = numerix.where(abs(P) < eps, eps, P)
        alpha = numerix.where(P > largeValue, (P - 1) / P, 0.5)
        Pmin = numerix.where(P > largeValue + 1, largeValue + 1, P)
        alpha = numerix.where((abs(Pmin) > eps) & (Pmin <= largeValue), ((Pmin - 1) * numerix.exp(Pmin) + 1) / (Pmin * (numerix.exp(Pmin) - 1)), alpha)

        return alpha

class ExponentialConvectionTerm(_AsymmetricConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::

       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the exponential scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """

    def _alpha(self, P):
        return _ExponentialConvectionTermAlpha(P)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
