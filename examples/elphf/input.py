#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "input.py"
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
This example adds two more components to
``examples/elphf/input1DphaseBinary.py``
one of which is another substitutional species and the other represents
electrons and diffuses interterstitially.

Parameters from `2004/January/21/elphf0214`

We start by defining a 1D mesh

    >>> from fipy import PhysicalField as PF

    >>> RT = (PF("1 Nav*kB") * PF("298 K"))
    >>> molarVolume = PF("1.80000006366754e-05 m**3/mol")
    >>> Faraday = PF("1 Nav*e")

    >>> L = PF("3 nm")
    >>> nx = 1200
    >>> dx = L / nx
    >>> # nx = 200
    >>> # dx = PF("0.01 nm")
    >>> ## dx = PF("0.001 nm") * (1.001 - 1/cosh(arange(-10, 10, .01)))
    >>> # L = nx * dx
    >>> mesh = Grid1D(dx = dx, nx = nx)
    >>> # mesh = Grid1D(dx = dx)
    >>> # L = mesh.facesRight[0].center[0] - mesh.facesLeft[0].center[0]
    >>> # L = mesh.cellCenters[0,-1] - mesh.cellCenters[0,0]


We create the phase field

    >>> timeStep = PF("1e-12 s")

    >>> phase = CellVariable(mesh = mesh, name = 'xi', value = 1, hasOld = 1)
    >>> phase.mobility = PF("1 m**3/J/s") / (molarVolume / (RT * timeStep))
    >>> phase.gradientEnergy = PF("3.6e-11 J/m") / (mesh.scale**2 * RT / molarVolume)

    >>> def p(xi):
    ...     return xi**3 * (6 * xi**2 - 15 * xi + 10.)

    >>> def g(xi):
    ...     return (xi * (1 - xi))**2

    >>> def pPrime(xi):
    ...     return 30. * (xi * (1 - xi))**2

    >>> def gPrime(xi):
    ...     return 4 * xi * (1 - xi) * (0.5 - xi)

We create four components

    >>> class ComponentVariable(CellVariable):
    ...     def __init__(self, mesh, value = 0., name = '', standardPotential = 0., barrier = 0., diffusivity = None, valence = 0, equation = None, hasOld = 1):
    ...         self.standardPotential = standardPotential
    ...         self.barrier = barrier
    ...         self.diffusivity = diffusivity
    ...         self.valence = valence
    ...         self.equation = equation
    ...         CellVariable.__init__(self, mesh = mesh, value = value, name = name, hasOld = hasOld)
    ...
    ...     def copy(self):
    ...         return self.__class__(mesh = self.mesh, value = self.value,
    ...                               name = self.name,
    ...                               standardPotential = self.standardPotential,
    ...                               barrier = self.barrier,
    ...                               diffusivity = self.diffusivity,
    ...                               valence = self.valence,
    ...                               equation = self.equation,
    ...                               hasOld = 0)

the solvent

    >>> solvent = ComponentVariable(mesh = mesh, name = 'H2O', value = 1.)
    >>> CnStandardPotential = PF("34139.7265625 J/mol") / RT
    >>> CnBarrier = PF("3.6e5 J/mol") / RT
    >>> CnValence = 0

and two solute species

    >>> substitutionals = [
    ...     ComponentVariable(mesh = mesh, name = 'SO4',
    ...                       diffusivity = PF("1e-9 m**2/s") / (mesh.scale**2/timeStep),
    ...                       standardPotential = PF("24276.6640625 J/mol") / RT,
    ...                       barrier = CnBarrier,
    ...                       valence = -2,
    ...                       value = PF("0.000010414586295976 mol/l") * molarVolume),
    ...     ComponentVariable(mesh = mesh, name = 'Cu',
    ...                       diffusivity = PF("1e-9 m**2/s") / (mesh.scale**2/timeStep),
    ...                       standardPotential = PF("-7231.81396484375 J/mol") / RT,
    ...                       barrier = CnBarrier,
    ...                       valence = +2,
    ...                       value = PF("55.5553718417909 mol/l") * molarVolume)]

and one interstitial

    >>> interstitials = [
    ...     ComponentVariable(mesh = mesh, name = 'e-',
    ...                       diffusivity = PF("1e-9 m**2/s") / (mesh.scale**2/timeStep),
    ...                       standardPotential = PF("-33225.9453125 J/mol") / RT,
    ...                       barrier = 0.,
    ...                       valence = -1,
    ...                       value = PF("111.110723815414 mol/l") * molarVolume)]

    >>> for component in substitutionals:
    ...     solvent -= component
    ...     component.standardPotential -= CnStandardPotential
    ...     component.barrier -= CnBarrier
    ...     component.valence -= CnValence

Finally, we create the electrostatic potential field

    >>> potential = CellVariable(mesh = mesh, name = 'phi', value = 0.)

    >>> permittivity = PF("78.49 eps0") / (Faraday**2 * mesh.scale**2 / (RT * molarVolume))

    >>> permittivity = 1.
    >>> permitivityPrime = 0.

The thermodynamic parameters are chosen to give a solid phase rich in electrons
and the solvent and a liquid phase rich in the two substitutional species

.. warning: Addition and subtraction cause `solvent` to lose some crucial information
   so we only append it after the fact.

..

    >>> solvent.standardPotential = CnStandardPotential
    >>> solvent.barrier = CnBarrier
    >>> solvent.valence = CnValence

Once again, we start with a sharp phase boundary

    >>> x = mesh.cellCenters[0]
    >>> phase.setValue(x < L / 2)
    >>> interstitials[0].setValue("0.000111111503177394 mol/l" * molarVolume, where=x > L / 2)
    >>> substitutionals[0].setValue("0.249944439430068 mol/l" * molarVolume, where=x > L / 2)
    >>> substitutionals[1].setValue("0.249999982581341 mol/l" * molarVolume, where=x > L / 2)

We again create the phase equation as in ``examples.elphf.phase.input1D``

    >>> mesh.setScale(1)

    >>> phase.equation = TransientTerm(coeff = 1/phase.mobility) \
    ...     == DiffusionTerm(coeff = phase.gradientEnergy) \
    ...     - (permitivityPrime / 2.) * potential.grad.dot(potential.grad)

We linearize the source term in the same way as in `example.phase.simple.input1D`.

    >>> enthalpy = solvent.standardPotential
    >>> barrier = solvent.barrier
    >>> for component in substitutionals + interstitials:
    ...     enthalpy += component * component.standardPotential
    ...     barrier += component * component.barrier

    >>> mXi = -(30 * phase * (1 - phase) * enthalpy +  4 * (0.5 - phase) * barrier)
    >>> dmXidXi = (-60 * (0.5 - phase) * enthalpy + 4 * barrier)
    >>> S1 = dmXidXi * phase * (1 - phase) + mXi * (1 - 2 * phase)
    >>> S0 = mXi * phase * (1 - phase) - phase * S1

    >>> phase.equation -= S0 + ImplicitSourceTerm(coeff = S1)

and we create the diffustion equation for the solute as in
``examples.elphf.diffusion.input1D``

    >>> for Cj in substitutionals:
    ...     CkSum = ComponentVariable(mesh = mesh, value = 0.)
    ...     CkFaceSum = FaceVariable(mesh = mesh, value = 0.)
    ...     for Ck in [Ck for Ck in substitutionals if Ck is not Cj]:
    ...         CkSum += Ck
    ...         CkFaceSum += Ck.harmonicFaceValue
    ...
    ...     counterDiffusion = CkSum.faceGrad
    ...     # phaseTransformation = (pPrime(phase.harmonicFaceValue) * Cj.standardPotential
    ...     #         + gPrime(phase.harmonicFaceValue) * Cj.barrier) * phase.faceGrad
    ...     phaseTransformation = (pPrime(phase).harmonicFaceValue * Cj.standardPotential
    ...             + gPrime(phase).harmonicFaceValue * Cj.barrier) * phase.faceGrad
    ...     # phaseTransformation = (p(phase).faceGrad * Cj.standardPotential
    ...     #         + g(phase).faceGrad * Cj.barrier)
    ...     electromigration = Cj.valence * potential.faceGrad
    ...     convectionCoeff = counterDiffusion + \
    ...         solvent.harmonicFaceValue * (phaseTransformation + electromigration)
    ...     convectionCoeff *= (Cj.diffusivity / (1. - CkFaceSum))
    ...
    ...     Cj.equation = (TransientTerm()
    ...                    == DiffusionTerm(coeff=Cj.diffusivity)
    ...                    + PowerLawConvectionTerm(coeff=convectionCoeff))

    >>> for Cj in interstitials:
    ...     # phaseTransformation = (pPrime(phase.harmonicFaceValue) * Cj.standardPotential
    ...     #         + gPrime(phase.harmonicFaceValue) * Cj.barrier) * phase.faceGrad
    ...     phaseTransformation = (pPrime(phase).harmonicFaceValue * Cj.standardPotential
    ...             + gPrime(phase).harmonicFaceValue * Cj.barrier) * phase.faceGrad
    ...     # phaseTransformation = (p(phase).faceGrad * Cj.standardPotential
    ...     #         + g(phase).faceGrad * Cj.barrier)
    ...     electromigration = Cj.valence * potential.faceGrad
    ...     convectionCoeff = Cj.diffusivity * (1 + Cj.harmonicFaceValue) * \
    ...         (phaseTransformation + electromigration)
    ...
    ...     Cj.equation = (TransientTerm()
    ...                    == DiffusionTerm(coeff=Cj.diffusivity)
    ...                    + PowerLawConvectionTerm(coeff=convectionCoeff))

And Poisson's equation

    >>> charge = 0.
    >>> for Cj in interstitials + substitutionals:
    ...     charge += Cj * Cj.valence

    >>> potential.equation = DiffusionTerm(coeff = permittivity) + charge == 0

If running interactively, we create viewers to display the results

    >>> if __name__ == '__main__':
    ...     phaseViewer = Viewer(vars=phase, datamin=0, datamax=1)
    ...     concViewer = Viewer(vars=[solvent] + substitutionals + interstitials, ylog=True)
    ...     potentialViewer = Viewer(vars = potential)
    ...     phaseViewer.plot()
    ...     concViewer.plot()
    ...     potentialViewer.plot()
    ...     raw_input("Press a key to continue")

Again, this problem does not have an analytical solution, so after
iterating to equilibrium

    >>> solver = LinearLUSolver(tolerance = 1e-3)

    >>> potential.constrain(0., mesh.facesLeft)

    >>> phase.residual = CellVariable(mesh = mesh)
    >>> potential.residual = CellVariable(mesh = mesh)
    >>> for Cj in substitutionals + interstitials:
    ...     Cj.residual = CellVariable(mesh = mesh)
    >>> residualViewer = Viewer(vars = [phase.residual, potential.residual] + [Cj.residual for Cj in substitutionals + interstitials])

    >>> tsv = TSVViewer(vars = [phase, potential] + substitutionals + interstitials)

    >>> dt = substitutionals[0].diffusivity * 100
    >>> # dt = 1.
    >>> elapsed = 0.
    >>> maxError = 1e-1
    >>> SAFETY = 0.9
    >>> ERRCON = 1.89e-4
    >>> desiredTimestep = 1.
    >>> thisTimeStep = 0.
    >>> print "%3s: %20s | %20s | %20s | %20s" % ("i", "elapsed", "this", "next dt", "residual")
    >>> residual = 0.
    >>> for i in range(500): # iterate
    ...     if thisTimeStep == 0.:
    ...         tsv.plot(filename = "%s.tsv" % str(elapsed * timeStep))
    ...
    ...     for field in [phase, potential] + substitutionals + interstitials:
    ...         field.updateOld()
    ...
    ...     while 1:
    ...         for j in range(10): # sweep
    ...             print i, j, dt * timeStep, residual
    ...             # raw_input()
    ...             residual = 0.
    ...
    ...             phase.equation.solve(var = phase, dt = dt)
    ...             # print phase.name, phase.equation.residual.max()
    ...             residual = max(phase.equation.residual.max(), residual)
    ...             phase.residual[:] = phase.equation.residual
    ...
    ...             potential.equation.solve(var = potential, dt = dt)
    ...             # print potential.name, potential.equation.residual.max()
    ...             residual = max(potential.equation.residual.max(), residual)
    ...             potential.residual[:] = potential.equation.residual
    ...
    ...             for Cj in substitutionals + interstitials:
    ...                 Cj.equation.solve(var = Cj,
    ...                                   dt = dt,
    ...                                   solver = solver)
    ...                 # print Cj.name, Cj.equation.residual.max()
    ...                 residual = max(Cj.equation.residual.max(), residual)
    ...                 Cj.residual[:] = Cj.equation.residual
    ...
    ...             # print
    ...             # phaseViewer.plot()
    ...             # concViewer.plot()
    ...             # potentialViewer.plot()
    ...             # residualViewer.plot()
    ...
    ...         residual /= maxError
    ...         if residual <= 1.:
    ...             break	# step succeeded
    ...
    ...         dt = max(SAFETY * dt * residual**-0.2, 0.1 * dt)
    ...         if thisTimeStep + dt == thisTimeStep:
    ...             raise FloatingPointError("step size underflow")
    ...
    ...     thisTimeStep += dt
    ...
    ...     if residual > ERRCON:
    ...         dt *= SAFETY * residual**-0.2
    ...     else:
    ...         dt *= 5.
    ...
    ...     # dt *= (maxError / residual)**0.5
    ...
    ...     if thisTimeStep >= desiredTimestep:
    ...         elapsed += thisTimeStep
    ...         thisTimeStep = 0.
    ...     else:
    ...         dt = min(dt, desiredTimestep - thisTimeStep)
    ...
    ...     if __name__ == '__main__':
    ...         phaseViewer.plot()
    ...         concViewer.plot()
    ...         potentialViewer.plot()
    ...         print "%3d: %20s | %20s | %20s | %g" % (i, str(elapsed * timeStep), str(thisTimeStep * timeStep), str(dt * timeStep), residual)


we confirm that the far-field phases have remained separated

    >>> ends = take(phase, (0,-1))
    >>> allclose(ends, (1.0, 0.0), rtol = 1e-5, atol = 1e-5)
    1

and that the concentration fields has appropriately segregated into into
their respective phases

    >>> ends = take(interstitials[0], (0,-1))
    >>> allclose(ends, (0.4, 0.3), rtol = 3e-3, atol = 3e-3)
    1
    >>> ends = take(substitutionals[0], (0,-1))
    >>> allclose(ends, (0.3, 0.4), rtol = 3e-3, atol = 3e-3)
    1
    >>> ends = take(substitutionals[1], (0,-1))
    >>> allclose(ends, (0.1, 0.2), rtol = 3e-3, atol = 3e-3)
    1
"""
__docformat__ = 'restructuredtext'

## def _test():
##     import doctest
##     return doctest.testmod()
##
## if __name__ == "__main__":
##     _test()
##     raw_input("finished")

if __name__ == '__main__':
    ## from fipy.tools.profiler.profiler import Profiler
    ## from fipy.tools.profiler.profiler import calibrate_profiler

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    # profile.stop()

    raw_input("finished")
