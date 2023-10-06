r"""
A simple 1D example to test the setup of the multicomponent diffusion
equations.  The diffusion equation for each species in single-phase
multicomponent system can be expressed as

.. math::

   \frac{\partial C_j}{\partial t}
   = D_{jj}\nabla^2 C_j
     + D_{j}\nabla\cdot
       \frac{C_j}{1 - \sum_{\substack{k=2\\ k \neq j}}^{n-1} C_k}
           \sum_{\substack{i=2\\ i \neq j}}^{n-1} \nabla C_i


where :math:`C_j` is the concentration of the :math:`j^\text{th}` species,
:math:`t` is time, :math:`D_{jj}` is the self-diffusion coefficient of the
:math:`j^\text{th}` species, and :math:`\sum_{\substack{i=2\\ i \neq j}}^{n-1}`
represents the summation over all substitutional species in the system,
excluding the solvent and the component of interest.

We solve the problem on a 1D mesh

>>> nx = 400
>>> dx = 0.01
>>> L = nx * dx

>>> from fipy import CellVariable, FaceVariable, Grid1D, TransientTerm, DiffusionTerm, PowerLawConvectionTerm, DefaultAsymmetricSolver, Viewer
>>> mesh = Grid1D(dx = dx, nx = nx)

One component in this ternary system will be designated the "solvent"

>>> class ComponentVariable(CellVariable):
...     def __init__(self, mesh, value = 0., name = '',
...                  standardPotential = 0., barrier = 0.,
...                  diffusivity = None, valence = 0, equation = None):
...         CellVariable.__init__(self, mesh = mesh, value = value,
...                               name = name)
...         self.standardPotential = standardPotential
...         self.barrier = barrier
...         self.diffusivity = diffusivity
...         self.valence = valence
...         self.equation = equation
... 
...     def copy(self):
...         return self.__class__(mesh = self.mesh,
...                               value = self.value,
...                               name = self.name,
...                               standardPotential =
...                                   self.standardPotential,
...                               barrier = self.barrier,
...                               diffusivity = self.diffusivity,
...                               valence = self.valence,
...                               equation = self.equation)

>>> solvent = ComponentVariable(mesh = mesh, name = 'Cn', value = 1.)

We can create an arbitrary number of components,
simply by providing a :class:`tuple` or :class:`list` of components

>>> substitutionals = [
...     ComponentVariable(mesh = mesh, name = 'C1', diffusivity = 1.,
...                       standardPotential = 1., barrier = 1.),
...     ComponentVariable(mesh = mesh, name = 'C2', diffusivity = 1.,
...                       standardPotential = 1., barrier = 1.),
...     ]

>>> interstitials = []

>>> for component in substitutionals:
...     solvent -= component

We separate the solution domain into two different concentration regimes

>>> x = mesh.cellCenters[0]
>>> substitutionals[0].setValue(0.3)
>>> substitutionals[0].setValue(0.6, where=x > L / 2)
>>> substitutionals[1].setValue(0.6)
>>> substitutionals[1].setValue(0.3, where=x > L / 2)

We create one diffusion equation for each substitutional component

>>> for Cj in substitutionals:
...     CkSum = ComponentVariable(mesh = mesh, value = 0.)
...     CkFaceSum = FaceVariable(mesh = mesh, value = 0.)
...     for Ck in [Ck for Ck in substitutionals if Ck is not Cj]:
...         CkSum += Ck
...         CkFaceSum += Ck.harmonicFaceValue
... 
...     convectionCoeff = CkSum.faceGrad \
...                       * (Cj.diffusivity / (1. - CkFaceSum))
... 
...     Cj.equation = (TransientTerm()
...                    == DiffusionTerm(coeff=Cj.diffusivity)
...                    + PowerLawConvectionTerm(coeff=convectionCoeff))
...     Cj.solver = DefaultAsymmetricSolver(precon=None, iterations=3200)

If we are running interactively, we create a viewer to see the results

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=[solvent] + substitutionals,
...                     datamin=0, datamax=1)
...     viewer.plot()

Now, we iterate the problem to equilibrium, plotting as we go

>>> from builtins import range
>>> for i in range(40):
...     for Cj in substitutionals:
...         Cj.equation.solve(var=Cj,
...                           dt=10000.,
...                           solver=Cj.solver)
...     if __name__ == '__main__':
...         viewer.plot()

Since there is nothing to maintain the concentration separation in this problem,
we verify that the concentrations have become uniform

.. note::

   Between `petsc=3.13.2=h82b89f7_0` and `petsc=3.13.4=h82b89f7_0`, PETSc
   ceased achieving 1e-7 tolerance when solving on 2 processors on Linux.
   Solving on macOS is OK. Solving on 1, 3, or 4 processors is OK.

>>> print(substitutionals[0].allclose(0.45, rtol = 2e-7, atol = 2e-7))
True
>>> print(substitutionals[1].allclose(0.45, rtol = 2e-7, atol = 2e-7))
True
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    ## from fipy.tools.profiler.profiler import Profiler
    ## from fipy.tools.profiler.profiler import calibrate_profiler

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    # profile.stop()

    input("finished")
