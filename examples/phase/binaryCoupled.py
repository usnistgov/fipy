"""Simultaneously solve a phase-field evolution and solute diffusion problem in one-dimension.

It is straightforward to extend a phase field model to include binary alloys.
As in :mod:`examples.phase.simple`, we will examine a 1D problem

.. index::
   single: Grid1D

>>> from fipy import CellVariable, Variable, Grid1D, TransientTerm, DiffusionTerm, ImplicitSourceTerm, LinearLUSolver, Viewer, DefaultAsymmetricSolver
>>> from fipy.tools import numerix

>>> nx = 400
>>> dx = 5e-6 # cm
>>> L = nx * dx
>>> mesh = Grid1D(dx=dx, nx=nx)

The Helmholtz free energy functional can be written as the integral
:cite:`BoettingerReview:2002` :cite:`McFaddenReview:2002` :cite:`Wheeler:1992`

.. math::

   \\mathcal{F}\\left(\\phi, C, T\\right)
   = \\int_\\mathcal{V} \\left\\{
       f(\\phi, C, T)
       + \\frac{\\kappa_\\phi}{2}\\abs{\\nabla\\phi}^2
       + \\frac{\\kappa_C}{2}\\abs{\\nabla C}^2
   \\right\\} dV

over the volume :math:`\\mathcal{V}` as a function of phase :math:`\\phi`
[#phiCoupled]_

.. index::
   single: CellVariable

>>> phase = CellVariable(name="phase", mesh=mesh, hasOld=1)

composition :math:`C`

>>> C = CellVariable(name="composition", mesh=mesh, hasOld=1)

and temperature :math:`T` [#TCoupled]_

.. index::
   single: Variable

>>> T = Variable(name="temperature")

Frequently, the gradient energy term in concentration is ignored and we
can derive governing equations

.. math::
   :label: eq:phase:binaryCoupled:phase

    \\frac{\\partial\\phi}{\\partial t}
    = M_\\phi \\left( \\kappa_\\phi \\nabla^2 \\phi
                   - \\frac{\\partial f}{\\partial \\phi} \\right)

for phase and

.. math::
   :label: eq:phase:binaryCoupled:diffusion

   \\frac{\\partial C}{\\partial t}
   = \\nabla\\cdot\\left( M_C \\nabla \\frac{\\partial f}{\\partial C} \\right)

for solute.

The free energy density :math:`f(\\phi, C, T)` can be constructed in many
different ways. One approach is to construct free energy densities for
each of the pure components, as functions of phase, *e.g.*

.. math::

   f_A(\\phi, T) = p(\\phi) f_A^S(T)
   + \\left(1 - p(\\phi)\\right) f_A^L(T) + \\frac{W_A}{2} g(\\phi)

where :math:`f_A^L(T)`, :math:`f_B^L(T)`, :math:`f_A^S(T)`, and :math:`f_B^S(T)`
are the free energy densities of the pure components. There are a
variety of choices for the interpolation function :math:`p(\\phi)` and the
barrier function :math:`g(\\phi)`,

such as those shown in :mod:`examples.phase.simple`

>>> def p(phi):
...     return phi**3 * (6 * phi**2 - 15 * phi + 10)

>>> def g(phi):
...     return (phi * (1 - phi))**2

The desired thermodynamic model can then be applied to obtain
:math:`f(\\phi, C, T)`, such as for a regular solution,

.. math::

   f(\\phi, C, T) &= (1 - C) f_A(\\phi, T) + C f_B(\\phi, T) \\\\
   &\\qquad + R T \\left[
       (1 - C) \\ln (1 - C) + C \\ln C
   \\right]
   + C (1 - C) \\left[
       \\Omega_S p(\\phi)
       + \\Omega_L \\left( 1 - p(\\phi) \\right)
   \\right]

where

>>> R = 8.314 # J / (mol K)

is the gas constant and :math:`\\Omega_S` and :math:`\\Omega_L` are the
regular solution interaction parameters for solid and liquid.

Another approach is useful when the free energy densities :math:`f^L(C, T)`
and :math:`f^S(C,T)` of the alloy in the solid and liquid phases are
known. This might be the case when the two different phases have
different thermodynamic models or when one or both is obtained from a
Calphad code. In this case, we can construct

.. math::

   f(\\phi, C, T) = p(\\phi) f^S(C,T)
   + \\left(1 - p(\\phi)\\right) f^L(C, T)
   + \\left[
       (1-C) \\frac{W_A}{2} + C \\frac{W_B}{2}
   \\right] g(\\phi).

When the thermodynamic models are the same in both phases, both
approaches should yield the same result.

We choose the first approach and make the simplifying assumptions of an
ideal solution and that

.. math::

   f_A^L(T) & = 0 \\\\
   f_A^S(T) - f_A^L(T) &= \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}

and likewise for component :math:`B`.

>>> LA = 2350. # J / cm**3
>>> LB = 1728. # J / cm**3
>>> TmA = 1728. # K
>>> TmB = 1358. # K

>>> enthalpyA = LA * (T - TmA) / TmA
>>> enthalpyB = LB * (T - TmB) / TmB

This relates the difference between the free energy densities of the
pure solid and pure liquid phases to the latent heat :math:`L_A` and the
pure component melting point :math:`T_M^A`, such that

.. math::

   f_A(\\phi, T) = \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A} p(\\phi)
   + \\frac{W_A}{2} g(\\phi).

With these assumptions

.. math::
   :label: eq:phase:binaryCoupled:phaseTransformation

   \\frac{\\partial f}{\\partial \\phi}
   &= (1-C) \\frac{\\partial f_A}{\\partial \\phi}
   + C \\frac{\\partial f_B}{\\partial \\phi} \\nonumber \\\\
   &= \\left\\{
       (1-C) \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}
       + C \\frac{L_B\\left(T - T_M^B\\right)}{T_M^B}
   \\right\\} p'(\\phi)
   + \\left\\{
     (1-C) \\frac{W_A}{2} + C \\frac{W_B}{2}
   \\right\\} g'(\\phi)

and

.. math::
   :label: eq:phase:binaryCoupled:chemicalPotential

   \\frac{\\partial f}{\\partial C}
   &= \\left[f_B(\\phi, T) + \\frac{R T}{V_m} \\ln C\\right]
   - \\left[f_A(\\phi, T) + \\frac{R T}{V_m} \\ln (1-C) \\right] \\nonumber \\\\
   &= \\left[\\mu_B(\\phi, C, T) - \\mu_A(\\phi, C, T) \\right] / V_m

where :math:`\\mu_A` and :math:`\\mu_B` are the classical chemical potentials
for the binary species. :math:`p'(\\phi)` and :math:`g'(\\phi)` are the
partial derivatives of of :math:`p` and :math:`g` with respect to :math:`\\phi`

>>> def pPrime(phi):
...     return 30. * g(phi)

>>> def gPrime(phi):
...     return 2. * phi * (1 - phi) * (1 - 2 * phi)

:math:`V_m` is the molar volume, which we take to be independent of
concentration and phase

>>> Vm = 7.42 # cm**3 / mol

On comparison with :mod:`examples.phase.simple`, we can see that the
present form of the phase field equation is identical to the one found
earlier, with the source now composed of the concentration-weighted average
of the source for either pure component. We let the pure component barriers
equal the previous value

>>> deltaA = deltaB = 1.5 * dx
>>> sigmaA = 3.7e-5 # J / cm**2
>>> sigmaB = 2.9e-5 # J / cm**2
>>> betaA = 0.33 # cm / (K s)
>>> betaB = 0.39 # cm / (K s)
>>> kappaA = 6 * sigmaA * deltaA # J / cm
>>> kappaB = 6 * sigmaB * deltaB # J / cm
>>> WA = 6 * sigmaA / deltaA # J / cm**3
>>> WB = 6 * sigmaB / deltaB # J / cm**3

and define the averages

>>> W = (1 - C) * WA / 2. + C * WB / 2.
>>> enthalpy = (1 - C) * enthalpyA + C * enthalpyB

We can now linearize the source exactly as before

>>> mPhi = -((1 - 2 * phase) * W + 30 * phase * (1 - phase) * enthalpy)
>>> dmPhidPhi = 2 * W - 30 * (1 - 2 * phase) * enthalpy
>>> S1 = dmPhidPhi * phase * (1 - phase) + mPhi * (1 - 2 * phase)
>>> S0 = mPhi * phase * (1 - phase) - S1 * phase

Using the same gradient energy coefficient and phase field mobility

>>> kappa = (1 - C) * kappaA + C * kappaB
>>> Mphi = TmA * betaA / (6 * LA * deltaA)

we define the phase field equation

>>> phaseEq = (TransientTerm(1/Mphi, var=phase) == DiffusionTerm(coeff=kappa, var=phase)
...            + S0 + ImplicitSourceTerm(coeff=S1, var=phase))

----

When coding explicitly, it is typical to simply write a function to
evaluate the chemical potentials :math:`\\mu_A` and :math:`\\mu_B` and then
perform the finite differences necessary to calculate their gradient and
divergence, e.g.,::

    def deltaChemPot(phase, C, T):
        return ((Vm * (enthalpyB * p(phase) + WA * g(phase)) + R * T * log(1 - C)) -
                (Vm * (enthalpyA * p(phase) + WA * g(phase)) + R * T * log(C)))

    for j in range(faces):
        flux[j] = ((Mc[j+.5] + Mc[j-.5]) / 2) \\
          * (deltaChemPot(phase[j+.5], C[j+.5], T) \\
            - deltaChemPot(phase[j-.5], C[j-.5], T)) / dx

    for j in range(cells):
        diffusion = (flux[j+.5] - flux[j-.5]) / dx

where we neglect the details of the outer boundaries (``j = 0`` and ``j = N``)
or exactly how to translate ``j+.5`` or ``j-.5`` into an array index,
much less the complexities of higher dimensions. FiPy can handle all of
these issues automatically, so we could just write::

    chemPotA = Vm * (enthalpyA * p(phase) + WA * g(phase)) + R * T * log(C)
    chemPotB = Vm * (enthalpyB * p(phase) + WB * g(phase)) + R * T * log(1-C)
    flux = Mc * (chemPotB - chemPotA).faceGrad
    eq = TransientTerm() == flux.divergence

Although the second syntax would essentially work as written, such an
explicit implementation would be very slow. In order to take advantage
of :term:`FiPy`'s implicit solvers, it is necessary to reduce
Eq. :eq:`eq:phase:binaryCoupled:diffusion` to the canonical form of
Eq. :eq:`num:gen`, hence we must expand
Eq. :eq:`eq:phase:binaryCoupled:chemicalPotential` as

.. math::

   \\frac{\\partial f}{\\partial C}
   = \\left[
       \\frac{L_B\\left(T - T_M^B\\right)}{T_M^B}
       - \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}
   \\right] p(\\phi)
   + \\frac{R T}{V_m} \\left[\\ln C - \\ln (1-C)\\right]
   + \\frac{W_B - W_A}{2} g(\\phi)

In either bulk phase, :math:`\\nabla p(\\phi) = \\nabla g(\\phi) = 0`, so
we can then reduce Eq. :eq:`eq:phase:binaryCoupled:diffusion` to

.. math::
   :label: eq:phase:binaryCoupled:diffusion:bulk

   \\frac{\\partial C}{\\partial t}
   &= \\nabla\\cdot\\left( M_C \\nabla \\left\\{
       \\frac{R T}{V_m} \\left[\\ln C - \\ln (1-C)\\right]
   \\right\\}
   \\right) \\nonumber \\\\
   &= \\nabla\\cdot\\left[
       \\frac{M_C R T}{C (1-C) V_m} \\nabla C
   \\right]

and, by comparison with Fick's second law

.. math::

   \\frac{\\partial C}{\\partial t}
   = \\nabla\\cdot\\left[D \\nabla C\\right],

we can associate the mobility :math:`M_C` with the intrinsic diffusivity :math:`D_C` by
:math:`M_C \\equiv D_C C (1-C) V_m / R T` and write Eq. :eq:`eq:phase:binaryCoupled:diffusion` as

.. math::
   :label: eq:phase:binaryCoupled:diffusion:canonical

   \\frac{\\partial C}{\\partial t}
   &= \\nabla\\cdot\\left( D_C \\nabla C \\right) \\nonumber \\\\
   &\\qquad + \\nabla\\cdot\\left(
   \\frac{D_C C (1 - C) V_m}{R T}
   \\left\\{
       \\left[
           \\frac{L_B\\left(T - T_M^B\\right)}{T_M^B}
           - \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}
       \\right] \\nabla p(\\phi)
       + \\frac{W_B - W_A}{2} \\nabla g(\\phi)
   \\right\\}
   \\right). \\\\
   &= \\nabla\\cdot\\left( D_C \\nabla C \\right) \\nonumber \\\\
   &\\qquad + \\nabla\\cdot\\left(
   \\frac{D_C C (1 - C) V_m}{R T}
   \\left\\{
       \\left[
           \\frac{L_B\\left(T - T_M^B\\right)}{T_M^B}
           - \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}
       \\right] p'(\\phi)
       + \\frac{W_B - W_A}{2} g'(\\phi)
   \\right\\} \\nabla \\phi
   \\right).

The first term is clearly a :class:`~fipy.terms.diffusionTerm.DiffusionTerm` in :math:`C`. The second is a
:class:`~fipy.terms.diffusionTerm.DiffusionTerm` in :math:`\\phi` with a diffusion coefficient

.. math::

   D_{\\phi}(C, \\phi) =
   \\frac{D_C C (1 - C) V_m}{R T}
   \\left\\{
       \\left[
           \\frac{L_B\\left(T - T_M^B\\right)}{T_M^B}
           - \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}
       \\right] p'(\\phi)
       + \\frac{W_B - W_A}{2} g'(\\phi)
   \\right\\},

such that

.. math::

   \\frac{\\partial C}{\\partial t}
   = \\nabla\\cdot\\left( D_C \\nabla C \\right) + \\nabla\\cdot\\left(D_\\phi \\nabla \\phi \\right)


or

>>> Dl = Variable(value=1e-5) # cm**2 / s
>>> Ds = Variable(value=1e-9) # cm**2 / s
>>> Dc = (Ds - Dl) * phase.arithmeticFaceValue + Dl

>>> Dphi = ((Dc * C.harmonicFaceValue * (1 - C.harmonicFaceValue) * Vm / (R * T))
...         * ((enthalpyB - enthalpyA) * pPrime(phase.arithmeticFaceValue)
...            + 0.5 * (WB - WA) * gPrime(phase.arithmeticFaceValue)))

>>> diffusionEq = (TransientTerm(var=C)
...                == DiffusionTerm(coeff=Dc, var=C)
...                + DiffusionTerm(coeff=Dphi, var=phase))

>>> eq = phaseEq & diffusionEq

----

We initialize the phase field to a step function in the middle of the domain

>>> phase.setValue(1.)
>>> phase.setValue(0., where=mesh.cellCenters[0] > L/2.)

and start with a uniform composition field :math:`C = 1/2`

>>> C.setValue(0.5)

In equilibrium, :math:`\\mu_A(0, C_L, T) = \\mu_A(1, C_S, T)` and
:math:`\\mu_B(0, C_L, T) = \\mu_B(1, C_S, T)` and, for ideal solutions, we can
deduce the liquidus and solidus compositions as

.. math::

   C_L &= \\frac{1 - \\exp\\left(-\\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}\\frac{V_m}{R T}\\right)}
   {\\exp\\left(-\\frac{L_B\\left(T - T_M^B\\right)}{T_M^B}\\frac{V_m}{R T}\\right)
   - \\exp\\left(-\\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}\\frac{V_m}{R T}\\right)} \\\\
   C_S &= \\exp\\left(-\\frac{L_B\\left(T - T_M^B\\right)}{T_M^B}\\frac{V_m}{R T}\\right) C_L

.. index::
   single: exp

>>> Cl = ((1. - numerix.exp(-enthalpyA * Vm / (R * T)))
...       / (numerix.exp(-enthalpyB * Vm / (R * T)) - numerix.exp(-enthalpyA * Vm / (R * T))))
>>> Cs = numerix.exp(-enthalpyB * Vm / (R * T)) * Cl

The phase fraction is predicted by the lever rule

>>> Cavg = C.cellVolumeAverage
>>> fraction = (Cl - Cavg) / (Cl - Cs)

For the special case of ``fraction = Cavg = 0.5``, a little bit of algebra
reveals that the temperature that leaves the phase fraction unchanged is
given by

>>> T.setValue((LA + LB) * TmA * TmB / (LA * TmB + LB * TmA))

In this simple, binary, ideal solution case, we can derive explicit
expressions for the solidus and liquidus compositions. In general, this may
not be possible or practical. In that event, the root-finding facilities in
SciPy can be used.

We'll need a function to return the two conditions for equilibrium

.. math::

   0 = \\mu_A(1, C_S, T) - \\mu_A(0, C_L, T) &=
   \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A} V_m
   + R T \\ln (1 - C_S) - R T \\ln (1 - C_L) \\\\
   0 = \\mu_B(1, C_S, T) - \\mu_B(0, C_L, T) &=
   \\frac{L_B\\left(T - T_M^B\\right)}{T_M^B} V_m
   + R T \\ln C_S - R T \\ln C_L

.. index::
   single: log
   single: array

>>> def equilibrium(C):
...     return [numerix.array(enthalpyA * Vm
...                           + R * T * numerix.log(1 - C[0])
...                           - R * T * numerix.log(1 - C[1])),
...             numerix.array(enthalpyB * Vm
...                           + R * T * numerix.log(C[0])
...                           - R * T * numerix.log(C[1]))]

and we'll have much better luck if we also supply the Jacobian

.. math::

   \\left[\\begin{matrix}
       \\frac{\\partial(\\mu_A^S - \\mu_A^L)}{\\partial C_S}
       & \\frac{\\partial(\\mu_A^S - \\mu_A^L)}{\\partial C_L} \\\\
       \\frac{\\partial(\\mu_B^S - \\mu_B^L)}{\\partial C_S}
       & \\frac{\\partial(\\mu_B^S - \\mu_B^L)}{\\partial C_L}
   \\end{matrix}\\right]
   =
   R T\\left[\\begin{matrix}
       -\\frac{1}{1-C_S} & \\frac{1}{1-C_L} \\\\
       \\frac{1}{C_S} & -\\frac{1}{C_L}
   \\end{matrix}\\right]

>>> def equilibriumJacobian(C):
...     return R * T * numerix.array([[-1. / (1 - C[0]), 1. / (1 - C[1])],
...                                   [ 1. / C[0],      -1. / C[1]]])

.. index::
   pair: module; scipy

>>> try:
...     from scipy.optimize import fsolve # doctest: +SCIPY
...     CsRoot, ClRoot = fsolve(func=equilibrium, x0=[0.5, 0.5],
...                             fprime=equilibriumJacobian) # doctest: +SCIPY
... except ImportError:
...     ClRoot = CsRoot = 0
...     print("The SciPy library is not available to calculate the solidus and \
... liquidus concentrations")

>>> print(Cl.allclose(ClRoot)) # doctest: +SCIPY
1
>>> print(Cs.allclose(CsRoot)) # doctest: +SCIPY
1

We plot the result against the sharp interface solution

>>> sharp = CellVariable(name="sharp", mesh=mesh)
>>> x = mesh.cellCenters[0]
>>> sharp.setValue(Cs, where=x < L * fraction)
>>> sharp.setValue(Cl, where=x >= L * fraction)

>>> elapsed = Variable(value=0.) # s

.. index::
   pair: module; fipy.viewers

>>> if __name__ == '__main__':
...     try:
...         from examples.phase.phaseViewer import PhaseViewer
...
...         viewer = PhaseViewer(phase=phase, C=C, sharp=sharp,
...                              elapsed=elapsed,
...                              L=L, deltaA=deltaA,
...                              tmin=1e-5, tmax=300 * 3600,
...                              datamin=0., datamax=1.)
...     except ImportError:
...         viewer = Viewer(vars=(phase, C, sharp),
...                         datamin=0., datamax=1.)
...     viewer.plot()

Because the phase field interface will not move, and because we've seen in
earlier examples that the diffusion problem is unconditionally stable, we
need take only one very large timestep to reach equilibrium

>>> dt = 1.e5

Because the phase field equation is coupled to the composition through
``enthalpy`` and ``W`` and the diffusion equation is coupled to the phase
field through ``phaseTransformationVelocity``, it is necessary sweep this
non-linear problem to convergence. We use the "residual" of the equations
(a measure of how well they think they have solved the given set of linear
equations) as a test for how long to sweep.

.. index::
   single: LinearLUSolver
   single: solve
   single: sweep

We now use the ":meth:`~fipy.terms.term.Term.sweep`" method instead of
":meth:`~fipy.terms.term.Term.solve`" because we require the residual.

>>> import fipy.solvers.solver
>>> from fipy.tools import parallelComm
>>> if parallelComm.Nproc > 1:
...     if fipy.solvers.solver_suite == 'petsc':
...         solver = DefaultAsymmetricSolver(tolerance=1e-10, precon='hypre')
...     elif fipy.solvers.solver_suite in ['trilinos', 'no-pysparse']:
...         # Trilinos scales by initial residual
...         # b-vector L2norm is ~1e15
...         solver = DefaultAsymmetricSolver(tolerance=1e-24)
... else:
...     solver = LinearLUSolver(tolerance=1e-10)

>>> phase.updateOld()
>>> C.updateOld()
>>> res = 1.
>>> initialRes = None
>>> sweep = 0

>>> while res > 1e-8 and sweep < 100:
...     res = eq.sweep(dt=dt, solver=solver)
...     if initialRes is None:
...         initialRes = res
...     res = res / initialRes
...     sweep += 1

>>> from fipy import input
>>> if __name__ == '__main__':
...     viewer.plot()
...     input("Stationary phase field. Press <return> to proceed...")

.. image:: /figures/examples/phase/binary/stationary.*
   :width: 90%
   :align: center
   :alt: phase and composition fields in equilibrium, compared with phase diagram concentrations

We verify that the bulk phases have shifted to the predicted solidus and
liquidus compositions

>>> X = mesh.faceCenters[0]
>>> print(Cs.allclose(C.faceValue[X.value==0], atol=1e-2))
True
>>> print(Cl.allclose(C.faceValue[X.value==L], atol=1e-2))
True

and that the phase fraction remains unchanged

>>> print(fraction.allclose(phase.cellVolumeAverage, atol=2e-4))
1

while conserving mass overall

>>> print(Cavg.allclose(0.5, atol=1e-8))
1

----

We now quench by ten degrees

>>> T.setValue(T() - 10.) # K

>>> sharp.setValue(Cs, where=x < L * fraction)
>>> sharp.setValue(Cl, where=x >= L * fraction)

Because this lower temperature will induce the phase interface to move
(solidify), we will need to take much smaller timesteps (the time scales of
diffusion and of phase transformation compete with each other).

The `CFL limit`_ requires that no interface should advect more than one grid
spacing in a timestep. We can get a rough idea for the maximum timestep we can
take by looking at the velocity of convection induced by phase transformation in
Eq. :eq:`eq:phase:binaryCoupled:diffusion:canonical` (even though there is no explicit
convection in the coupled form used for this example, the principle remains the
same). If we assume that the phase changes from 1 to 0 in a single grid spacing,
that the diffusivity is `Dl` at the interface, and that the term due to the difference in
barrier heights is negligible:

.. math::

   \\vec{u}_\\phi &= \\frac{D_\\phi}{C} \\nabla \\phi
   \\\\
   &\\approx
   \\frac{D_l \\frac{1}{2} V_m}{R T}
   \\left[
       \\frac{L_B\\left(T - T_M^B\\right)}{T_M^B}
       - \\frac{L_A\\left(T - T_M^A\\right)}{T_M^A}
   \\right] \\frac{1}{\\Delta x}
   \\\\
   &\\approx
   \\frac{D_l \\frac{1}{2} V_m}{R T}
   \\left(L_B + L_A\\right) \\frac{T_M^A - T_M^B}{T_M^A + T_M^B}
   \\frac{1}{\\Delta x}
   \\\\
   &\\approx 0.28~\\mathrm{cm/s}

To get a :math:`\\text{CFL} = \\vec{u}_\\phi \\Delta t / \\Delta x < 1`, we need a
time step of about :math:`10^{-5}~\\mathrm{s}`.

>>> dt0 = 1.e-5

>>> from builtins import range
>>> for i in range(8):
...     phase.updateOld()
...     C.updateOld()
...     res = 1e+10
...     sweep = 0
...     while (res > 1e-3 or abs(Cavg.value - 0.5) > 1e-8) and sweep < 20:
...         res = eq.sweep(dt=dt0, solver=solver)
...         sweep += 1
...     elapsed.value = (i + 1) * dt0
...     if __name__ == '__main__':
...         viewer.plot()

>>> from fipy import input
>>> if __name__ == '__main__':
...     input("Moving phase field. Press <return> to proceed...")

.. image:: /figures/examples/phase/binary/binaryCoupled-0.000080.*
   :width: 90%
   :align: center

We see that the composition on either side of the interface approaches the
sharp-interface solidus and liquidus, but it will take a great many more
timesteps to reach equilibrium. If we waited sufficiently long, we
could again verify the final concentrations and phase fraction against the
expected values.

We can estimate the time to equilibration by examining the time for the
diffusion field to become uniform.  In the liquid, this will take
:math:`\\mathcal{O}((10~\\mathrm{\\mu m})^2 / D_l) =
0.1~\\mathrm{s}` and in the solid
:math:`\\mathcal{O}((10~\\mathrm{\\mu m})^2 / D_s) =
1000~\\mathrm{s}`.

Not wanting to take a hundred-million steps, we employ adaptive time
stepping, using the :term:`steppyingstounes` package.  This package takes
care of many of the messy details of stepping, like overshoot, underflow,
and step size adaptation, while keeping the structure of our solve loop
largely intact.

>>> from steppyngstounes import SequenceStepper, PIDStepper # doctest: +STEPPYNGSTOUNES
>>> from itertools import count

Assuming the process is dominated by diffusion, we can take steps that
increase geometrically.  Since we're unsure if diffusion is the only
process controlling dynamics, we take each increasing step with an adaptive
stepper that uses a `PID controller`_ to keep the equation residuals and
mass conservation within acceptable limits.  The total number of solves is
not strongly sensitive to the number of sweeps, but two sweeps seems to be
both sufficient and efficient.

We'll only advance the step if it's successful, so we need to update the
old values before we get started.

>>> phase.updateOld()
>>> C.updateOld()

>>> if __name__ == '__main__':
...     totaltime = 300 * 3600 # 300 h
... else:
...     totaltime = 32e-5 # 320 us

>>> dt = dt0

>>> mass_tolerance = 1e-6
>>> residual_tolerance = 1e-3

>>> if ((parallelComm.Nproc > 1)
...     and (fipy.solvers.solver_suite in ['trilinos', 'no-pysparse'])):
...     # Trilinos on linux in parallel doesn't conserve as well
...     mass_tolerance = 1e-5

>>> for checkpoint in SequenceStepper(start=float(elapsed), stop=totaltime,     # doctest: +STEPPYNGSTOUNES
...                                   sizes=(dt0 * 2**(n/2) for n in count(7))):
...     for step in PIDStepper(start=checkpoint.begin,
...                            stop=checkpoint.end,
...                            size=dt):
...         for sweep in range(2):
...             res = eq.sweep(dt=step.size, solver=solver)
...             # print(step.begin, step.size, sweep, res)
...         err = max(res / residual_tolerance,
...                   abs(Cavg.value - 0.5) / mass_tolerance)
...         if step.succeeded(error=err):
...             phase.updateOld()
...             C.updateOld()
...             elapsed.value = step.end
...             if __name__ == '__main__':
...                 viewer.plot()
...         else:
...             phase.value = phase.old
...             C.value = C.old
...     # the last step might have been smaller than possible,
...     # if it was near the end of the checkpoint range
...     dt = step.want
...     _ = checkpoint.succeeded()
...     if __name__ == '__main__':
...         viewer.plot()

>>> from fipy import input
>>> if __name__ == '__main__':
...     input("Re-equilbrated phase field. Press <return> to proceed...")

.. image:: /figures/examples/phase/binary/binaryCoupled-0.000899.*
   :width: 30%
   :alt: phase and composition fields at t=0.000899, compared with final phase diagram concentrations

.. image:: /figures/examples/phase/binary/binaryCoupled-8.949963.*
   :width: 30%
   :alt: phase and composition fields at t=8.949963, compared with final phase diagram concentrations

.. image:: /figures/examples/phase/binary/binaryCoupled-1080000.000000.*
   :width: 30%
   :alt: phase and composition fields at t=1080000, compared with final phase diagram concentrations

The interface moves :math:`\\approx 2.8~\\mathrm{\\mu m}` in
:math:`70~\\mathrm{ms}`, driven by diffusion in the liquid
phase (compare the estimate above of :math:`0.1~\\mathrm{s}`).
For the next
:math:`12~\\mathrm{s}`, the interface stalls while the solute step
trapped in the solid phase diffuses outward
(:math:`(2.8~\\mathrm{\\mu m})^2 / D_s =
\mathcal{O}(80~\\mathrm{s})`).  Once the solute gradient in the
solid reaches the new position of the interface, the solidification front
begins to move, driven by diffusion in the solid.  When the solute in the
solid becomes uniform, the interface stalls again after :math:`\\approx
4000~\\mathrm{s}`, having moved another
:math:`2.9~\\mathrm{\\mu m}` (recall the estimate of
:math:`1000~\\mathrm{s}` for equilibration in the solid).  After this
point, there is essentially no further motion of the interface and barely
perceptible changes in the concentration field.

.. note::

   This evolution is qualitatively consistent with that seen in
   :mod:`examples.phase.binary`, but the interface does not move as far and
   the bulk concentrations are further from the phase diagram values.  The
   computation also takes substantially longer than the uncoupled variant.

.. rubric:: Footnotes

.. [#phiCoupled] We will find that we need to "sweep" this non-linear problem
   (see *e.g.* the composition-dependent diffusivity example in
   :mod:`examples.diffusion.mesh1D`), so we declare :math:`\\phi` and :math:`C`
   to retain an "old" value.

.. [#TCoupled] we are going to want to
   examine different temperatures in this example, so we declare :math:`T`
   as a :class:`~fipy.variables.variable.Variable`

.. _CFL limit: http://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition

.. _PID controller: https://en.wikipedia.org/wiki/PID_controller
"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
