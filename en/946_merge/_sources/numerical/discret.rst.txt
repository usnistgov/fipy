Finite Volume Method
--------------------

To use the FVM, the solution domain must first be divided into
non-overlapping polyhedral elements or cells.  A solution domain
divided in such a way is generally known as a mesh (as we will see, a
:class:`~fipy.meshes.mesh.Mesh` is also a :term:`FiPy` object).  A mesh consists of vertices,
faces and cells (see Figure :ref:`fig:meshcartoon`).  In the FVM the
variables of interest are averaged over control volumes (CVs).  The
CVs are either defined by the cells or are centered on the vertices.

.. _fig:meshcartoon:

.. figure:: meshcartoon.*
   :width: 50%
   :align: center

   Mesh

   A mesh consists of cells, faces and vertices. For the purposes of
   :term:`FiPy`, the divider between two cells is known as a face for all
   dimensions.

Cell Centered FVM (CC-FVM)
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the CC-FVM the CVs are formed by the mesh cells with the cell
center "storing" the average variable value in the CV, (see
Figure :ref:`fig:vc-cc-fv`). The face fluxes are approximated using the
variable values in the two adjacent cells surrounding the face. This
low order approximation has the advantage of being efficient and
requiring matrices of low band width (the band width is equal to the
number of cell neighbors plus one) and thus low storage
requirement. However, the mesh topology is restricted due to
orthogonality and conjunctionality requirements. The value at a face is
assumed to be the average value over the face. On an unstructured mesh
the face center may not lie on the line joining the CV centers, which
will lead to an error in the face interpolation. :term:`FiPy` currently
only uses the CC-FVM.

Boundary Conditions
'''''''''''''''''''

The natural boundary condition for CC-FVM is no-flux.  For :eq:`num:gen`,
the boundary condition is

.. math::

   \hat{n}\cdot[\vec{u}\phi - (\Gamma_i\nabla)^n] = 0

Vertex Centered FVM (VC-FVM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the VC-FVM, the CV is centered around the vertices and the cells
are divided into sub-control volumes that make up the main CVs (see
Figure :ref:`fig:vc-cc-fv`). The vertices "store" the average
variable values over the CVs. The CV faces are constructed within the
cells rather than using the cell faces as in the CC-FVM. The face
fluxes use all the vertex values from the cell where the face is
located to calculate interpolations. For this reason, the VC-FVM is
less efficient and requires more storage (a larger matrix band width)
than the CC-FVM. However, the mesh topology does not have the same
restrictions as the CC-FVM. :term:`FiPy` does not have a VC-FVM capability.

.. Future releases of :term:`FiPy` will have both
   the CC-FVM and VC-FVM capabilities.

.. _fig:vc-cc-fv:

.. figure:: vc-cc-fv.*
   :width: 50%
   :align: center

   CV structure for an unstructured mesh

   (a) :math:`\Omega_a` represents a vertex-based CV and (b)
   :math:`\Omega_1`, :math:`\Omega_2`, :math:`\Omega_3` and
   :math:`\Omega_4` represent cell centered CVs.

.. _section:discretization:

Discretization
--------------

The first step in the discretization of Equation :eq:`num:gen`
using the CC-FVM is to integrate over a CV and then make appropriate
approximations for fluxes across the boundary of each CV.  In this
section, each term in Equation :eq:`num:gen` will be examined
separately.

Transient Term :math:`\partial (\rho \phi) / \partial t`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the transient term, the discretization of the integral :math:`\int_V` over the volume of a CV is given by

.. math::
   :label: num:tra

   \int_V \frac{\partial (\rho \phi)}{\partial t} dV
   \simeq
   \frac{(\rho_{P} \phi_{P} - \rho_{P}^\text{old} \phi_{P}^\text{old}) V_P}{\Delta t}

where :math:`\phi_P` represents the average value of :math:`\phi` in a CV
centered on a point :math:`P` and the superscript ":math:`\text{old}`"
represents the previous time-step value.  The value :math:`V_P` is the
volume of the CV and :math:`\Delta t` is the time step size.

This term is represented in :term:`FiPy` as

>>> TransientTerm(coeff=rho)

Convection Term :math:`\nabla \cdot \left( \vec{u} \phi \right)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The discretization for the convection term is given by

.. math::
   :label: num:con

   \int_V \nabla \cdot (\vec{u} \phi)\,dV &=
   \int_S (\vec{n} \cdot \vec{u})\phi\,dS \\
   &\simeq \sum_{f} (\vec{n} \cdot \vec{u})_f \phi_f A_f

where we have used the divergence theorem to transform the integral
over the CV volume :math:`\int_V` into an integral over the CV surface
:math:`\int_S`.  The summation over the faces of a CV is denoted by
:math:`\sum_{f}` and :math:`A_f` is the area of each face.  The vector :math:`\vec{n}`
is the normal to the face pointing out of the CV into an adjacent CV
centered on point :math:`A`.  When using a first order approximation,
the value of :math:`\phi_f` must depend on the average value in adjacent
cell :math:`\phi_A` and the average value in the cell of interest :math:`\phi_P`,
such that

.. math::

   \phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A.

The weighting factor :math:`\alpha_f` is determined by the convection
scheme, described in :ref:`sec:NumericalSchemes`.

This term is represented in :term:`FiPy` as

.. currentmodule:: fipy.terms

>>> <SpecificConvectionTerm>(coeff=u)

where :samp:`{<SpecificConvectionTerm>}` can be any of
:class:`~centralDiffConvectionTerm.CentralDifferenceConvectionTerm`,
:class:`~exponentialConvectionTerm.ExponentialConvectionTerm`,
:class:`~hybridConvectionTerm.HybridConvectionTerm`,
:class:`~powerLawConvectionTerm.PowerLawConvectionTerm`,
:class:`~upwindConvectionTerm.UpwindConvectionTerm`,
:class:`~explicitUpwindConvectionTerm.ExplicitUpwindConvectionTerm`, or
:class:`~vanLeerConvectionTerm.VanLeerConvectionTerm`.
The differences between these convection schemes are described
in Section :ref:`sec:NumericalSchemes`. The velocity coefficient
``u`` must be a rank-1 :class:`~fipy.variables.faceVariable.FaceVariable`, or a
constant vector in the form of a :term:`Python` list or
tuple, *e.g.* ``((1,), (2,))`` for a vector in 2D.

Diffusion Term :math:`\nabla \cdot \left( \Gamma_1 \nabla \phi \right)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The discretization for the diffusion term is given by

.. math::
   :label: num:dif

   \int_V \nabla \cdot (\Gamma\nabla\{\ldots\}) dV
   &= \int_S \Gamma (\vec{n} \cdot \nabla\{\ldots\}) dS \\
   &\simeq \sum_f \Gamma_f (\vec{n} \cdot \nabla\{\ldots\})_f A_f

:math:`\{\ldots\}` indicates recursive application of the specified
operation on :math:`\phi`, depending on
the order of the diffusion term.
The estimation for the flux, :math:`(\vec{n} \cdot \nabla\{\ldots\})_f`, is
obtained via

.. math::

   (\vec{n} \cdot \nabla\{\ldots\})_f \simeq \frac{\{\ldots\}_A-\{\ldots\}_P}{d_{AP}}

where the value of :math:`d_{AP}` is the distance between neighboring cell
centers.  This estimate relies on the orthogonality of the mesh, and
becomes increasingly inaccurate as the non-orthogonality increases.
Correction terms have been derived to improve this error but are not
currently included in :term:`FiPy` :cite:`croftphd`.

This term is represented in :term:`FiPy` as

.. currentmodule:: fipy.terms

>>> DiffusionTerm(coeff=Gamma1)

or

>>> ExplicitDiffusionTerm(coeff=Gamma1)

:class:`~explicitDiffusionTerm.ExplicitDiffusionTerm` is provided primarily for
illustrative purposes, although :mod:`examples.diffusion.mesh1D`
demonstrates its use in Crank-Nicolson time stepping.
:class:`~implicitDiffusionTerm.ImplicitDiffusionTerm` is almost always
preferred (:class:`~diffusionTerm.DiffusionTerm` is a synonym for
:class:`~implicitDiffusionTerm.ImplicitDiffusionTerm` to reinforce this
preference). One can also create an explicit diffusion term with

>>> (Gamma1 * phi.faceGrad).divergence

.. _discret-higherOrderDiffusion:

Higher Order Diffusion
''''''''''''''''''''''

Higher order diffusion expressions, such as :math:`\nabla^4 \phi` or
:math:`\nabla \cdot \left( \Gamma_1 \nabla \left( \nabla\cdot\left(
\Gamma_2 \nabla \phi\right) \right) \right)` for  Cahn-Hilliard are
represented as

>>> DiffusionTerm(coeff=(Gamma1, Gamma2))

The number of elements supplied for ``coeff`` determines the
order of the term.

.. note::

   While this multiple-coefficient form is still supported,
   :ref:`CoupledEquations` are the recommended approach for higher order
   expressions.

Source Term
~~~~~~~~~~~

Any term that cannot be written in one of the previous forms is considered
a source :math:`S_{\phi}`. The discretization for the source term is given
by,

.. math::
   :label: num:sou

   \int_V S_{\phi}\,dV \simeq S_\phi V_P.

Including any negative dependence of :math:`S_\phi` on :math:`\phi` increases
solution stability. The dependence can only be included in a linear
manner so Equation :eq:`num:sou` becomes

.. math::

   V_P (S_0 + S_1 \phi_P),

where :math:`S_0` is the source which is independent of :math:`\phi` and
:math:`S_1` is the coefficient of the source which is linearly dependent
on :math:`\phi`.

A source term is represented in :term:`FiPy` essentially as it appears in
mathematical form, *e.g.*, :math:`3\kappa^2 + b \sin
\theta` would be written

>>> 3 * kappa**2 + b * numerix.sin(theta)

.. note::

   Functions like :func:`~numpy.sin` can be obtained from the
   :mod:`fipy.tools.numerix` module.

   .. warning::

      Generally, things will not work as expected if the equivalent
      function is used from the :term:`NumPy` or :term:`SciPy` library.

If, however, the source depends on the variable that is being solved for,
it can be advantageous to linearize the source and cast part of it as an
implicit source term, *e.g.*, :math:`3\kappa^2 + \phi \sin \theta`
might be written as

>>> 3 * kappa**2 + ImplicitSourceTerm(coeff=sin(theta))

.. _section:linear-equations:

Linear Equations
----------------

The aim of the discretization is to reduce the continuous general
equation to a set of discrete linear equations that can then be solved
to obtain the value of the dependent variable at each CV center. This
results in a sparse linear system that requires an efficient iterative
scheme to solve. The iterative schemes available to :term:`FiPy` are
currently encapsulated in the :term:`Pysparse` and :term:`PyTrilinos`
suites of solvers and include most common solvers such as the conjugate
gradient method and LU decomposition.

Combining Equations :eq:`num:tra`, :eq:`num:con`,
:eq:`num:dif` and :eq:`num:sou`, the complete
discretization for equation :eq:`num:gen` can now be written for
each CV as

.. math::
   :label: num:dis

   \frac{\rho_{P}(\phi_{P} - \phi_P^\text{old}) V_P}{\Delta t}
   + &
   \sum_{f} (\vec{n} \cdot \vec{u})_f A_f
   \left[\alpha_f \phi_P +\left(1-\alpha_f\right)\phi_A\right]
   \nonumber \\
   = &
   \sum_f \Gamma_f A_f \frac{(\phi_A-\phi_P)}{d_{AP}}
   +
   V_P ( S_0 + S_1 \phi_P ).


Equation :eq:`num:dis` is now in the form of a set of linear
combinations between each CV value and its neighboring values and can be
written in the form

.. math::
   :label: num:dap

   a_P \phi_P = \sum_f a_{A} \phi_{A} + b_P,

where

.. math::

   a_P &=  \frac{\rho_P V_P}{\Delta t} + \sum_f (a_{A} +
   F_f) - V_P S_1, \\
   a_{A} &=  D_f - ( 1 - \alpha_f ) F_f, \\
   b_P &= V_P S_0 + \frac{\rho_P V_P \phi_P^\text{old}}{\Delta t}.

The face coefficients, :math:`F_f` and :math:`D_f`, represent the convective strength
and diffusive conductance respectively, and are given by

.. math::

   F_f &= A_f ( \vec{u} \cdot \vec{n} )_f, \\
   D_f &= \frac{A_f \Gamma_f}{d_{AP}}.
