#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "surfactantBulkDiffusionEquation.py"
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

from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
from fipy.variables.levelSetDiffusionVariable import _LevelSetDiffusionVariable
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.diffusionTerm import DiffusionTermNoCorrection

def buildSurfactantBulkDiffusionEquation(bulkVar = None,
                                         distanceVar = None,
                                         surfactantVar = None,
                                         otherSurfactantVar = None,
                                         diffusionCoeff = None,
                                         transientCoeff = 1.,
                                         rateConstant = None):
    r"""

    The `buildSurfactantBulkDiffusionEquation` function returns a bulk diffusion of a
    species with a source term for the jump from the bulk to an interface.
    The governing equation is given by,

    .. math::

       \frac{\partial c}{\partial t} = \nabla \cdot D \nabla  c

    where,

    .. math::

       D = \begin{cases}
           D_c & \text{when $\phi > 0$} \\
           0  & \text{when $\phi \le 0$}
       \end{cases}

    The jump condition at the interface is defined by Langmuir
    adsorption. Langmuir adsorption essentially states that the ability for
    a species to jump from an electrolyte to an interface is proportional to
    the concentration in the electrolyte, available site density and a
    jump coefficient. The boundary condition at the interface is given by

    .. math::

       D \hat{n} \cdot \nabla c = -k c (1 - \theta) \qquad \text{at $\phi = 0$}.

    :Parameters:
      - `bulkVar`: The bulk surfactant concentration variable.
      - `distanceVar`: A `DistanceVariable` object
      - `surfactantVar`: A `SurfactantVariable` object
      - `otherSurfactantVar`: Any other surfactants that may remove this one.
      - `diffusionCoeff`: A float or a `FaceVariable`.
      - `transientCoeff`: In general 1 is used.
      - `rateConstant`: The adsorption coefficient.

    """

    spCoeff = rateConstant * distanceVar.cellInterfaceAreas / bulkVar.mesh.cellVolumes
    spSourceTerm = ImplicitSourceTerm(spCoeff)

    bulkSpCoeff = spCoeff * bulkVar
    coeff = bulkSpCoeff * surfactantVar.interfaceVar

    diffusionCoeff = _LevelSetDiffusionVariable(distanceVar,
                                                diffusionCoeff)

    eq =  TransientTerm(transientCoeff) - DiffusionTermNoCorrection(diffusionCoeff)

    if otherSurfactantVar is not None:
        otherCoeff = bulkSpCoeff * otherSurfactantVar.interfaceVar
    else:
        otherCoeff = 0

    return eq - coeff + spSourceTerm - otherCoeff
