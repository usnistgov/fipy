#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "explicitUpwindConvectionTerm.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.terms.abstractUpwindConvectionTerm import _AbstractUpwindConvectionTerm
from fipy.tools import numerix
from fipy.terms import TransientTermError

__all__ = ["ExplicitUpwindConvectionTerm"]

class ExplicitUpwindConvectionTerm(_AbstractUpwindConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::

       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P^\text{old} +(1-\alpha_f)\phi_A^\text{old}` and
    :math:`\alpha_f` is calculated using the upwind scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """

    def _getOldAdjacentValues(self, oldArray, id1, id2, dt):
        if dt is None:
            raise TransientTermError
        return numerix.take(oldArray, id1), numerix.take(oldArray, id2)

    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        weight = _AbstractUpwindConvectionTerm._getWeight(self, var, transientGeomCoeff, diffusionGeomCoeff)
        if 'implicit' in weight.keys():
            weight['explicit'] = weight['implicit']
            del weight['implicit']

        return weight
