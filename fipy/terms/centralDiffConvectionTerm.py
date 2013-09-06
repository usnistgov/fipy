#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "centralDiffConvectionTerm.py"
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.terms.abstractConvectionTerm import _AbstractConvectionTerm
from fipy.variables.faceVariable import FaceVariable

__all__ = ["CentralDifferenceConvectionTerm"]

class _CentralDifferenceConvectionTermAlpha(FaceVariable):
    def __init__(self, P):
        FaceVariable.__init__(self, P.mesh, value=0.5)

class CentralDifferenceConvectionTerm(_AbstractConvectionTerm):
    r"""

    This :class:`~fipy.terms.term.Term` represents

    .. math::
    
       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f
       
    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the central differencing scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """
    def _alpha(self, P):
        return _CentralDifferenceConvectionTermAlpha(P)

    
