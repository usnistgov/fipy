from __future__ import division
from __future__ import unicode_literals
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

    Parameters
    ----------
    bulkVar : ~fipy.variables.cellVariable.CellVariable
        The bulk surfactant concentration variable.
    distanceVar : ~fipy.variables.distanceVariable.DistanceVariable
    surfactantVar : ~fipy.variables.surfactantVariable.SurfactantVariable
    otherSurfactantVar : ~fipy.variables.surfactantVariable.SurfactantVariable
        Any other surfactants that may remove this one.
    diffusionCoeff : float or ~fipy.variables.faceVariable.FaceVariable
    transientCoeff : float
        In general 1 is used.
    rateConstant : float
        The adsorption coefficient.

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
