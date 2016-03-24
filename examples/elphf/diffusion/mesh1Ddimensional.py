#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mesh1Ddimensional.py"
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
 #
 # ###################################################################
 ##

r"""
In this example, we present the same three-component diffusion problem
introduced in ``examples/elphf/diffusion/mesh1D.py``
but we demonstrate FiPy's facility to use dimensional quantities.

    >>> import warnings
    >>> warnings.warn("\n\n\tSupport for physical dimensions is incomplete.\n\tIt is not possible to solve dimensional equations.\n")

    >>> from fipy import CellVariable, FaceVariable, PhysicalField, Grid1D, TransientTerm, DiffusionTerm, PowerLawConvectionTerm, LinearLUSolver, Viewer
    >>> from fipy.tools import numerix

We solve the problem on a 40 mm long 1D mesh

    >>> nx = 40
    >>> dx = PhysicalField(1.,"mm")
    >>> L = nx * dx
    >>> mesh = Grid1D(dx = dx, nx = nx)

Again, one component in this ternary system will be designated the "solvent"

    >>> class ComponentVariable(CellVariable):
    ...     def __init__(self, mesh, value = 0., name = '',
    ...                  standardPotential = 0., barrier = 0.,
    ...                  diffusivity = None, valence = 0, equation = None):
    ...         CellVariable.__init__(self, mesh = mesh, value = value,
    ...                               name = name)
    ...         self.standardPotential = Variable(standardPotential)
    ...         self.barrier = Variable(barrier)
    ...         self.diffusivity = Variable(diffusivity)
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

    >>> solvent = ComponentVariable(mesh = mesh, name = 'Cn', value = "1 mol/m**3")

We can create an arbitrary number of components,
simply by providing a `Tuple` or `list` of components

    >>> substitutionals = [
    ...     ComponentVariable(mesh = mesh, name = 'C1', diffusivity = "1e-9 m**2/s",
    ...                       standardPotential = 1., barrier = 1., value = "0.3 mol/m**3"),
    ...     ComponentVariable(mesh = mesh, name = 'C2', diffusivity = "1e-9 m**2/s",
    ...                       standardPotential = 1., barrier = 1., value = "0.6 mol/m**3"),
    ...     ]

    >>> interstitials = []

    >>> for component in substitutionals:
    ...     solvent -= component

We separate the solution domain into two different concentration regimes

    >>> x = mesh.cellCenters[0]
    >>> substitutionals[0].setValue("0.3 mol/m**3")
    >>> substitutionals[0].setValue("0.6 mol/m**3", where=x > L / 2)
    >>> substitutionals[1].setValue("0.6 mol/m**3")
    >>> substitutionals[1].setValue("0.3 mol/m**3", where=x > L / 2)

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
    ...                    + PowerLawConvectionTerm(coeff = convectionCoeff))

If we are running interactively, we create a viewer to see the results

    >>> if __name__ == '__main__':
    ...     viewer = Viewer(vars=[solvent] + substitutionals,
    ...                     datamin=0, datamax=1)
    ...     viewer.plot()

Now, we iterate the problem to equilibrium, plotting as we go

    >>> solver = LinearLUSolver()

    >>> for i in range(40):
    ...     for Cj in substitutionals:
    ...         Cj.updateOld()
    ...     for Cj in substitutionals:
    ...         Cj.equation.solve(var = Cj,
    ...                           dt = "1000 s",
    ...                           solver = solver)
    ...     if __name__ == '__main__':
    ...         viewer.plot()

Since there is nothing to maintain the concentration separation in this problem,
we verify that the concentrations have become uniform

    >>> print substitutionals[0].scaled.allclose("0.45 mol/m**3",
    ...     atol = "1e-7 mol/m**3", rtol = 1e-7)
    1
    >>> print substitutionals[1].scaled.allclose("0.45 mol/m**3",
    ...     atol = "1e-7 mol/m**3", rtol = 1e-7)
    1

.. note::

   The absolute tolerance `atol` must be in units compatible with the value to
   be checked, but the relative tolerance `rtol` is dimensionless.
"""
__docformat__ = 'restructuredtext'


if __name__ == '__main__':
    ## from fipy.tools.profiler.profiler import Profiler
    ## from fipy.tools.profiler.profiler import calibrate_profiler

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    # profile.stop()

    input("finished")
