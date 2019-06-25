from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'


from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
from fipy.variables.levelSetDiffusionVariable import _LevelSetDiffusionVariable
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.diffusionTerm import DiffusionTermNoCorrection

def buildMetalIonDiffusionEquation(ionVar = None,
                                   distanceVar = None,
                                   depositionRate = 1,
                                   transientCoeff = 1,
                                   diffusionCoeff = 1,
                                   metalIonMolarVolume = 1):

    r"""

    The `MetalIonDiffusionEquation` solves the diffusion of the metal
    species with a source term at the electrolyte interface. The governing
    equation is given by,

    .. math::

       \frac{\partial c}{\partial t} = \nabla \cdot D \nabla  c

    where,

    .. math::

       D = \begin{cases}
           D_c & \text{when $\phi > 0$} \\
           0  & \text{when $\phi \le 0$}
       \end{cases}

    The velocity of the interface generally has a linear dependence on ion
    concentration. The following boundary condition applies at the zero
    level set,

    .. math::

       D \hat{n} \cdot \nabla c = \frac{v(c)}{\Omega} \qquad \text{at $phi = 0$}

    where

    .. math::

       v(c) = c V_0

    The test case below is for a 1D steady state problem. The solution is
    given by:

    .. math::

       c(x) = \frac{c^{\infty}}{\Omega D / V_0 + L}\left(x - L\right) + c^{\infty}

    This is the test case,

    >>> from fipy.meshes import Grid1D
    >>> nx = 11
    >>> dx = 1.
    >>> from fipy.tools import serialComm
    >>> mesh = Grid1D(nx = nx, dx = dx, communicator=serialComm)
    >>> x, = mesh.cellCenters
    >>> from fipy.variables.cellVariable import CellVariable
    >>> ionVar = CellVariable(mesh = mesh, value = 1.)
    >>> from fipy.variables.distanceVariable \
    ...     import DistanceVariable
    >>> disVar = DistanceVariable(mesh = mesh,
    ...                           value = (x - 0.5) - 0.99,
    ...                           hasOld = 1)

    >>> v = 1.
    >>> diffusion = 1.
    >>> omega = 1.
    >>> cinf = 1.

    >>> eqn = buildMetalIonDiffusionEquation(ionVar = ionVar,
    ...                                      distanceVar = disVar,
    ...                                      depositionRate = v * ionVar,
    ...                                      diffusionCoeff = diffusion,
    ...                                      metalIonMolarVolume = omega)

    >>> ionVar.constrain(cinf, mesh.facesRight)

    >>> from builtins import range
    >>> for i in range(10):
    ...     eqn.solve(ionVar, dt = 1000)

    >>> L = (nx - 1) * dx - dx / 2
    >>> gradient = cinf / (omega * diffusion / v + L)
    >>> answer = gradient * (x - L - dx * 3 / 2) + cinf
    >>> answer[x < dx] = 1
    >>> print(ionVar.allclose(answer))
    1

    Testing the interface source term

    >>> from fipy.meshes import Grid2D
    >>> from fipy import numerix, serialComm
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2, communicator=serialComm)
    >>> from fipy.variables.distanceVariable import DistanceVariable
    >>> distance = DistanceVariable(mesh = mesh, value = (-.5, .5, .5, 1.5))
    >>> ionVar = CellVariable(mesh = mesh, value = (1, 1, 1, 1))
    >>> depositionRate = CellVariable(mesh=mesh, value=(1, 1, 1, 1))
    >>> source = depositionRate * distance.cellInterfaceAreas / mesh.cellVolumes / ionVar
    >>> sqrt = numerix.sqrt(2)
    >>> ans = CellVariable(mesh=mesh, value=(0, 1 / sqrt, 1 / sqrt, 0))
    >>> print(numerix.allclose(source, ans))
    True
    >>> distance[:] = (-1.5, -0.5, -0.5, 0.5)
    >>> print(numerix.allclose(source, (0, 0, 0, sqrt)))
    True

    Parameters
    ----------
    ionVar : ~fipy.variables.cellVariable.CellVariable
         The metal ion concentration variable.
    distanceVar : ~fipy.variables.distanceVariable.DistanceVariable
    depositionRate : float or ~fipy.variables.cellVariable.CellVariable
         Represents the interface deposition rate.
    transientCoeff : float
         The transient coefficient.
    diffusionCoeff : float or ~fipy.variables.faceVariable.FaceVariable
         The diffusion coefficient
    metalIonMolarVolume : float
         Molar volume of the metal ions.

    """

    diffusionCoeff = _LevelSetDiffusionVariable(distanceVar,
                                                diffusionCoeff)

    eq =  TransientTerm(transientCoeff) - DiffusionTermNoCorrection(diffusionCoeff)

    mesh = distanceVar.mesh
    coeff = depositionRate * distanceVar.cellInterfaceAreas / (mesh.cellVolumes * metalIonMolarVolume) / ionVar

    return eq + ImplicitSourceTerm(coeff)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
