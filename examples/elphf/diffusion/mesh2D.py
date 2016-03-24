#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mesh2D.py"
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
The same three-component diffusion problem as introduced in :mod:`examples/elphf/diffusion/mesh1D.py` but in 2D:
    >>> from fipy import CellVariable, FaceVariable, Grid2D, TransientTerm, DiffusionTerm, PowerLawConvectionTerm, Viewer

    >>> nx = 40
    >>> dx = 1.
    >>> L = nx * dx
    >>> mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = nx)

One component in this ternary system will be designated the "solvent"

    >>> class ComponentVariable(CellVariable):
    ...     def __init__(self, mesh, value = 0., name = '', standardPotential = 0.,
    ...                  barrier = 0., diffusivity = None, valence = 0, equation = None):
    ...         CellVariable.__init__(self, mesh = mesh, value = value, name = name)
    ...         self.standardPotential = standardPotential
    ...         self.barrier = barrier
    ...         self.diffusivity = diffusivity
    ...         self.valence = valence
    ...         self.equation = equation
    ...
    ...     def copy(self):
    ...         return self.__class__(mesh = self.mesh, value = self.value,
    ...                               name = self.name,
    ...                               standardPotential = self.standardPotential,
    ...                               barrier = self.barrier,
    ...                               diffusivity = self.diffusivity,
    ...                               valence = self.valence,
    ...                               equation = self.equation)

    >>> solvent = ComponentVariable(mesh = mesh, name = 'Cn', value = 1.)

We can create an arbitrary number of components,
simply by providing a `Tuple` or `list` of components

    >>> substitutionals = [
    ...     ComponentVariable(mesh = mesh, name = 'C1', diffusivity = 1.,
    ...                       standardPotential = 1., barrier = 1.),
    ...     ComponentVariable(mesh = mesh, name = 'C2', diffusivity = 1.,
    ...                       standardPotential = 1., barrier = 1.),
    ...     ]

    >>> interstitials = []

    >>> for component in substitutionals:
    ...     solvent -= component

Although we are not interested in them for this problem, we create one field to represent the "phase" (1 everywhere)

    >>> phase = CellVariable(mesh = mesh, name = 'xi', value = 1.)

and one field to represent the electrostatic potential (0 everywhere)

    >>> potential = CellVariable(mesh = mesh, name = 'phi', value = 0.)

Althought it is constant in this problem, in later problems we will need the following
functions of the phase field

    >>> def pPrime(xi):
    ...     return 30. * (xi * (1 - xi))**2

    >>> def gPrime(xi):
    ...     return 2 * xi * (1 - xi) * (1 - 2 * xi)

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
    ...     counterDiffusion = CkSum.faceGrad
    ...     phaseTransformation = \
    ...         (pPrime(phase.harmonicFaceValue) * Cj.standardPotential \
    ...         + gPrime(phase.harmonicFaceValue) * Cj.barrier) \
    ...             * phase.faceGrad
    ...     electromigration = Cj.valence * potential.faceGrad
    ...     convectionCoeff = counterDiffusion \
    ...         + solvent.harmonicFaceValue \
    ...             * (phaseTransformation + electromigration)
    ...     convectionCoeff *= (Cj.diffusivity / (1. - CkFaceSum))
    ...
    ...     Cj.equation = (TransientTerm()
    ...                    == DiffusionTerm(coeff=Cj.diffusivity)
    ...                    + PowerLawConvectionTerm(coeff=convectionCoeff))

If we are running interactively, we create a viewer to see the results

    >>> if __name__ == '__main__':
    ...     viewers = [Viewer(vars=field, datamin=0, datamax=1)
    ...                for field in [solvent] + substitutionals]
    ...     for viewer in viewers:
    ...         viewer.plot()
    ...     steps = 40
    ...     tol = 1e-7
    ... else:
    ...     steps = 20
    ...     tol = 1e-4

Now, we iterate the problem to equilibrium, plotting as we go

    >>> for i in range(steps):
    ...     for Cj in substitutionals:
    ...         Cj.equation.solve(var = Cj,
    ...                           dt = 10000)
    ...     if __name__ == '__main__':
    ...         for viewer in viewers:
    ...             viewer.plot()

Since there is nothing to maintain the concentration separation in this problem,
we verify that the concentrations have become uniform

    >>> substitutionals[0].allclose(0.45, rtol = tol, atol = tol).value
    1
    >>> substitutionals[1].allclose(0.45, rtol = tol, atol = tol).value
    1

We now rerun the problem with an initial condition that only has a
concentration step in one corner.

    >>> x, y = mesh.cellCenters
    >>> substitutionals[0].setValue(0.3)
    >>> substitutionals[0].setValue(0.6, where=(x > L / 2.) & (y > L / 2.))
    >>> substitutionals[1].setValue(0.6)
    >>> substitutionals[1].setValue(0.3, where=(x > L / 2.) & (y > L / 2.))

We iterate the problem to equilibrium again

    >>> for i in range(steps):
    ...     for Cj in substitutionals:
    ...         Cj.equation.solve(var = Cj,
    ...                           dt = 10000)
    ...     if __name__ == '__main__':
    ...         for viewer in viewers:
    ...             viewer.plot()

and verify that the correct uniform concentrations are achieved

    >>> substitutionals[0].allclose(0.375, rtol = tol, atol = tol).value
    1
    >>> substitutionals[1].allclose(0.525, rtol = tol, atol = tol).value
    1


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
