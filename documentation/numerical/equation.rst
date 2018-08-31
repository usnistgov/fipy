.. _sec:GeneralConservationEquation:

General Conservation Equation
-----------------------------

The equations that model the evolution of physical, chemical and
biological systems often have a remarkably universal form.  Indeed,
PDEs have proven necessary to model complex physical systems and
processes that involve variations in both space and time.  In general,
given a variable of interest :math:`\phi` such as species concentration,
pH, or temperature, there exists an evolution equation of the form

.. math::
   :label: general-equation

   \frac{\partial \phi}{\partial t} = H(\phi, \lambda_i)

where :math:`H` is a function of :math:`\phi`, other state variables
:math:`\lambda_i`, and higher order derivatives of all of these variables.
Examples of such systems are wide ranging, but include problems that
exhibit a combination of diffusing and reacting species, as well as
such diverse problems as determination of the electric potential in
heart tissue, of fluid flow, stress evolution, and even the
|Schroedinger| equation.

A general conservation equation, solved using :term:`FiPy`, can include any
combination of the following terms,

.. math::
   :label: num:gen

   \underbrace{
     \frac{\partial (\rho \phi)}{\partial t}
   }_{\text{transient}}
   +
   \underbrace{
     \vphantom{\frac{\partial (\rho \phi)}{\partial t}}
     \nabla \cdot \left( \vec{u} \phi \right)
   }_{\text{convection}}
   =
   \underbrace{
     \vphantom{\frac{\partial (\rho \phi)}{\partial t}}
     \left[ \nabla \cdot \left( \Gamma_i \nabla \right) \right]^n \phi
   }_{\text{diffusion}}
   +
   \underbrace{
     \vphantom{\frac{\partial (\rho \phi)}{\partial t}}
     S_{\phi}
   }_{\text{source}}

where :math:`\rho`, :math:`\vec{u}` and :math:`\Gamma_i` represent coefficients in the
transient, convection and diffusion terms, respectively.  These
coefficients can be arbitrary functions of any parameters or variables
in the system.  The variable :math:`\phi` represents the unknown quantity in
the equation.  The diffusion term can represent any higher order
diffusion-like term, where the order is given by the exponent :math:`n`.
For example, the diffusion term can represent conventional Fickian
diffusion [*i.e.*, :math:`\nabla\cdot(\Gamma\nabla\phi)`] when the
exponent :math:`n = 1` or a Cahn-Hilliard term [*i.e.*, :math:`\nabla
\cdot (\Gamma_1 \nabla [ \nabla \cdot \Gamma_2 \nabla \phi ) ] )`
:cite:`CahnHilliardI` :cite:`CahnHilliardII` :cite:`CahnHilliardIII`] when :math:`n = 2`,
or a phase field crystal term [*i.e.*, :math:`\nabla \cdot (\Gamma_1 \nabla
[ \nabla \cdot \Gamma_2 \nabla \{ \nabla \cdot \Gamma_3 \nabla \phi ) \} )
] )` :cite:`Elder:2011p2811`] when :math:`n = 3`, although spectral methods are probably a
better approach. Higher order terms (:math:`n > 3`) are also possible, but
the matrix condition number becomes quite poor.

.. |Schroedinger| unicode:: Schr U+00F6 dinger
