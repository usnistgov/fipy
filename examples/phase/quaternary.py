#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "quaternary.py"
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

r""" Solve a phase-field evolution and diffusion of four species in one-dimension.

The same procedure used to construct the two-component phase field
diffusion problem in :mod:`examples.phase.binary` can be used to build up a
system of multiple components. Once again, we'll focus on 1D.

.. index:: Grid1D

>>> from fipy import CellVariable, Grid1D, TransientTerm, DiffusionTerm, ImplicitSourceTerm, PowerLawConvectionTerm, DefaultAsymmetricSolver, Viewer
>>> from fipy.tools import numerix

>>> nx = 400
>>> dx = 0.01
>>> L = nx * dx
>>> mesh = Grid1D(dx = dx, nx = nx)

We consider a free energy density :math:`f(\phi, C_0,\ldots,C_N, T)`
that is a function of phase :math:`\phi`

.. index:: CellVariable

>>> phase = CellVariable(mesh=mesh, name='phase', value=1., hasOld=1)

interstitial components :math:`C_0 \ldots C_M`

>>> interstitials = [
...     CellVariable(mesh=mesh, name='C0', hasOld=1)
... ]

substitutional components :math:`C_{M+1} \ldots C_{N-1}`

>>> substitutionals = [
...     CellVariable(mesh=mesh, name='C1', hasOld=1),
...     CellVariable(mesh=mesh, name='C2', hasOld=1),
... ]

a "solvent" :math:`C_N` that is constrained by the concentrations of the
other substitutional species, such that :math:`C_N = 1 - \sum_{j=M}^{N-1}
C_j`,

>>> solvent = 1
>>> for Cj in substitutionals:
...     solvent -= Cj
>>> solvent.name = 'CN'

and temperature :math:`T`

>>> T = 1000

The free energy density of such a system can be written as

.. math::

   f(\phi, C_0, \ldots, C_N, T)
   &= \sum_{j=0}^N C_j \left[ \mu^\circ_j(\phi, T) + R T \ln \frac{C_j}{\rho} \right]

where

>>> R = 8.314 # J / (mol K)

is the gas constant. As in the binary case,

.. math::

   \mu^\circ_j(\phi, T) = p(\phi) \mu_j^{\circ S}(T)
   + \left(1 - p(\phi)\right) \mu_j^{\circ L}(T) + \frac{W_j}{2} g(\phi)

is constructed with the free energies of the pure components in each
phase, given the "tilting" function

>>> def p(phi):
...     return phi**3 * (6 * phi**2 - 15 * phi + 10)

and the "double well" function

>>> def g(phi):
...     return (phi * (1 - phi))**2

We consider a very simplified model that has partial molar volumes
:math:`\bar{V}_0 = \cdots = \bar{V}_{M} = 0` for the "interstitials" and
:math:`\bar{V}_{M+1} = \cdots = \bar{V}_{N} = 1` for the "substitutionals".
This approximation has been used in a number of models where density
effects are ignored, including the treatment of electrons in
electrodeposition processes :cite:`ElPhFI` :cite:`ElPhFII`. Under these
constraints

.. math::

   \frac{\partial f}{\partial \phi}
   &= \sum_{j=0}^N C_j \frac{\partial f_j}{\partial \phi}
   \nonumber \\
   &= \sum_{j=0}^N C_j \left[
       \mu_j^{\circ SL}(T) p'(\phi) + \frac{W_j}{2} g'(\phi)
   \right]
   \\
   \frac{\partial f}{\partial C_j}
   &= \left[\mu^\circ_j(\phi, T) + R T \ln \frac{C_j}{\rho} \right]
   \nonumber \\
   &= \mu_j(\phi, C_j , T)
   \qquad\text{for \( j = 0\ldots M \)}

and

.. math::

   \frac{\partial f}{\partial C_j}
   &= \left[\mu^\circ_j(\phi, T) + R T \ln \frac{C_j}{\rho} \right]
   - \left[\mu^\circ_N(\phi, T) + R T \ln \frac{C_N}{\rho} \right]
   \nonumber \\
   &= \left[\mu_j(\phi, C_j, T) - \mu_N(\phi, C_N, T) \right]
   \qquad\text{for \( j = M+1\ldots N-1 \)}

where :math:`\mu_j^{\circ SL}(T) \equiv \mu_j^{\circ S}(T) - \mu_j^{\circ
L}(T)` and where :math:`\mu_j` is the classical chemical potential of
component :math:`j` for the binary species and :math:`\rho = 1 +
\sum_{j=0}^{M} C_j` is the total molar density.

>>> rho = 1.
>>> for Cj in interstitials:
...     rho += Cj

:math:`p'(\phi)` and :math:`g'(\phi)` are the partial derivatives of of :math:`p`
and :math:`g` with respect to :math:`\phi`

>>> def pPrime(phi):
...     return 30. * g(phi)

>>> def gPrime(phi):
...     return 2. * phi * (1 - phi) * (1 - 2 * phi)

We "cook" the standard potentials to give the desired solid and liquid
concentrations, with a solid phase rich in interstitials and the solvent
and a liquid phase rich in the two substitutional species.

>>> interstitials[0].S = 0.3
>>> interstitials[0].L = 0.4
>>> substitutionals[0].S = 0.4
>>> substitutionals[0].L = 0.3
>>> substitutionals[1].S = 0.2
>>> substitutionals[1].L = 0.1
>>> solvent.S = 1.
>>> solvent.L = 1.
>>> for Cj in substitutionals:
...     solvent.S -= Cj.S
...     solvent.L -= Cj.L

>>> rhoS = rhoL = 1.
>>> for Cj in interstitials:
...     rhoS += Cj.S
...     rhoL += Cj.L

.. index:: log

>>> for Cj in interstitials + substitutionals + [solvent]:
...     Cj.standardPotential = R * T * (numerix.log(Cj.L/rhoL)
...                                     - numerix.log(Cj.S/rhoS))

>>> for Cj in interstitials:
...     Cj.diffusivity = 1.
...     Cj.barrier = 0.

>>> for Cj in substitutionals:
...     Cj.diffusivity = 1.
...     Cj.barrier = R * T

>>> solvent.barrier = R * T

-----

We create the phase equation

.. math::

   \frac{1}{M_\phi}\frac{\partial \phi}{\partial t}
   = \kappa_\phi \nabla^2 \phi
   - \sum_{j=0}^N C_j \left[
              \mu_j^{\circ SL}(T) p'(\phi) + \frac{W_j}{2} g'(\phi)
          \right]

with a semi-implicit source just as in :mod:`examples.phase.simple` and
:mod:`examples.phase.binary`

>>> enthalpy = 0.
>>> barrier = 0.
>>> for Cj in interstitials + substitutionals + [solvent]:
...     enthalpy += Cj * Cj.standardPotential
...     barrier += Cj * Cj.barrier

>>> mPhi = -((1 - 2 * phase) * barrier + 30 * phase * (1 - phase) * enthalpy)
>>> dmPhidPhi = 2 * barrier - 30 * (1 - 2 * phase) * enthalpy
>>> S1 = dmPhidPhi * phase * (1 - phase) + mPhi * (1 - 2 * phase)
>>> S0 = mPhi * phase * (1 - phase) - S1 * phase

>>> phase.mobility = 1.
>>> phase.gradientEnergy = 25
>>> phase.equation = TransientTerm(coeff=1/phase.mobility) \
...   == DiffusionTerm(coeff=phase.gradientEnergy) \
...      + S0 + ImplicitSourceTerm(coeff = S1)

We could construct the diffusion equations one-by-one, in the manner of
:mod:`examples.phase.binary`, but it is better to take advantage of the full
scripting power of the Python language, where we can easily loop over
components or even make "factory" functions if we desire. For the
interstitial diffusion equations, we arrange in canonical form as before:

.. math::

   \underbrace{
       \frac{\partial C_j}{\partial t}
       \vphantom{\left\{
           \overbrace{
               \left[
                   \mu_j^{\circ SL} \nabla p(\phi)
               \right]
           }^{\text{phase transformation}}
       \right\}}
   }_{\text{transient}}
   &= \underbrace{
       D_j\nabla^2 C_j
       \vphantom{\left\{
           \overbrace{
               \left[
                   \mu_j^{\circ SL} \nabla p(\phi)
               \right]
           }^{\text{phase transformation}}
       \right\}}
   }_{\text{diffusion}} \\
   & \qquad + \underbrace{
       D_j\nabla\cdot
       \frac{C_j}{1 + \sum_{\substack{k=0\\ k \neq j}}^{M} C_k}
       \left\{
           \overbrace{
               \frac{\rho}{R T}
               \left[
                   \mu_j^{\circ SL} \nabla p(\phi)
                   + \frac{W_j}{2} \nabla g(\phi)
               \right]
           }^{\text{phase transformation}}
           -
           \overbrace{
               \sum_{\substack{i=0\\ i \neq j}}^{M} \nabla C_i
           }^{\text{counter diffusion}}
       \right\}
   }_{\text{convection}}


.. index:: PowerLawConvectionTerm

>>> for Cj in interstitials:
...     phaseTransformation = (rho.harmonicFaceValue / (R * T)) \
...       * (Cj.standardPotential * p(phase).faceGrad
...          + 0.5 * Cj.barrier * g(phase).faceGrad)
...
...     CkSum = CellVariable(mesh=mesh, value=0.)
...     for Ck in [Ck for Ck in interstitials if Ck is not Cj]:
...         CkSum += Ck
...
...     counterDiffusion = CkSum.faceGrad
...
...     convectionCoeff = counterDiffusion + phaseTransformation
...     convectionCoeff *= (Cj.diffusivity
...                         / (1. + CkSum.harmonicFaceValue))
...
...     Cj.equation = (TransientTerm()
...                    == DiffusionTerm(coeff=Cj.diffusivity)
...                    + PowerLawConvectionTerm(coeff=convectionCoeff))

-----

The canonical form of the substitutional diffusion equations is

.. math::

   \underbrace{
        \frac{\partial C_j}{\partial t}
   }_{\text{transient}}
    &= \underbrace{
        D_{j}\nabla^2 C_j
        \vphantom{\frac{\partial C_j}{\partial t}}
    }_{\text{diffusion}} \\
    & \qquad + \underbrace{
        D_{j}\nabla\cdot
        \frac{C_j}{1 - \sum_{\substack{k=M+1\\ k \neq j}}^{N-1} C_k}
        \left\{
           \overbrace{
                \frac{C_N}{R T}
                \left[
                    \left(\mu_j^{\circ SL} - \mu_N^{\circ SL}\right) \nabla p(\phi)
                    + \frac{W_j - W_N}{2} \nabla g(\phi)
                \right]
                \vphantom{\sum_{\substack{i=M+1\\ i \neq j}}^{N-1} \nabla C_i}
           }^{\text{phase transformation}}
           +
           \overbrace{
               \sum_{\substack{i=M+1\\ i \neq j}}^{N-1} \nabla C_i
           }^{\text{counter diffusion}}
        \right\}
    }_{\text{convection}}

>>> for Cj in substitutionals:
...     phaseTransformation = (solvent.harmonicFaceValue / (R * T)) \
...       * ((Cj.standardPotential - solvent.standardPotential) * p(phase).faceGrad
...          + 0.5 * (Cj.barrier - solvent.barrier) * g(phase).faceGrad)
...
...     CkSum = CellVariable(mesh=mesh, value=0.)
...     for Ck in [Ck for Ck in substitutionals if Ck is not Cj]:
...         CkSum += Ck
...
...     counterDiffusion = CkSum.faceGrad
...
...     convectionCoeff = counterDiffusion + phaseTransformation
...     convectionCoeff *= (Cj.diffusivity
...                         / (1. - CkSum.harmonicFaceValue))
...
...     Cj.equation = (TransientTerm()
...                    == DiffusionTerm(coeff=Cj.diffusivity)
...                    + PowerLawConvectionTerm(coeff=convectionCoeff))

-----

We start with a sharp phase boundary

.. math::

   \xi =
   \begin{cases}
       1& \text{for $x \le L/2$,} \\
       0& \text{for $x > L/2$,}
   \end{cases}

>>> x = mesh.cellCenters[0]
>>> phase.setValue(1.)
>>> phase.setValue(0., where=x > L / 2)

and with uniform concentration fields, initially equal to the average of
the solidus and liquidus concentrations

>>> for Cj in interstitials + substitutionals:
...     Cj.setValue((Cj.S + Cj.L) / 2.)

If we're running interactively, we create a viewer

.. index::
   :module: viewers

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=([phase]
...                           + interstitials + substitutionals
...                           + [solvent]),
...                     datamin=0, datamax=1)
...     viewer.plot()

and again iterate to equilibrium

.. .. index:: DefaultAsymmetricSolver

>>> solver = DefaultAsymmetricSolver(tolerance=1e-10)


>>> dt = 10000
>>> for i in range(5):
...     for field in [phase] + substitutionals + interstitials:
...         field.updateOld()
...     phase.equation.solve(var = phase, dt = dt)
...     for field in substitutionals + interstitials:
...         field.equation.solve(var = field,
...                              dt = dt,
...                              solver = solver)
...     if __name__ == '__main__':
...         viewer.plot()

.. image:: quaternary.*
   :width: 90%
   :align: center
   :alt: phase and four composition fields in equilibrium

We can confirm that the far-field phases have remained separated

.. .. index:: allclose

>>> X = mesh.faceCenters[0]
>>> print numerix.allclose(phase.faceValue[X.value==0], 1.0, rtol = 1e-5, atol = 1e-5)
True
>>> print numerix.allclose(phase.faceValue[X.value==L], 0.0, rtol = 1e-5, atol = 1e-5)
True

and that the concentration fields have appropriately segregated into
their equilibrium values in each phase

>>> equilibrium = True
>>> for Cj in interstitials + substitutionals:
...     equilibrium &= numerix.allclose(Cj.faceValue[X.value==0], Cj.S, rtol = 3e-3, atol = 3e-3).value
...     equilibrium &= numerix.allclose(Cj.faceValue[X.value==L], Cj.L, rtol = 3e-3, atol = 3e-3).value
>>> print equilibrium
True

.. .. bibmissing:: /documentation/refs.bib
    :sort:
"""
__docformat__ = 'restructuredtext'

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input("finished")

## if __name__ == '__main__':
##     ## from fipy.tools.profiler.profiler import Profiler
##     ## from fipy.tools.profiler.profiler import calibrate_profiler
##
##     # fudge = calibrate_profiler(10000)
##     # profile = Profiler('profile', fudge=fudge)
##
##     import fipy.tests.doctestPlus
##     exec(fipy.tests.doctestPlus.getScript())
##
##     # profile.stop()
##
##     raw_input("finished")
##
