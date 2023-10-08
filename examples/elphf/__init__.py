r"""The following examples exhibit various parts of a model to study
electrochemical interfaces.  In a pair of papers, Guyer, Boettinger, Warren
and McFadden :cite:`ElPhFI` :cite:`ElPhFII` have shown that an electrochemical
interface can be modeled by an equation for the phase field :math:`\xi`

.. math::
   :label: elphf:phase

    \underbrace{
	\frac{1}{M_\xi}\frac{\partial \xi}{\partial t}
	\vphantom{\sum_{j=1}^{n} C_j}
    }_{\text{transient}}
    = 
    \underbrace{
	\kappa_{\xi}\nabla^2 \xi
	\vphantom{\sum_{j=1}^{n} C_j}
    }_{\text{diffusion}}
    - 
    \underbrace{
	\overbrace{
	    \sum_{j=1}^{n} C_j \left[
		p'(\xi) \Delta\mu_j^\circ
		+ g'(\xi) W_j
	    \right]
	}^{\text{phase transformation}}
	+
	\overbrace{
	    \frac{\epsilon'(\xi)}{2}\left(\nabla\phi\right)^2
	    \vphantom{\sum_{j=1}^{n} C_j}
	}^{\text{dielectric}}
    }_{\text{source}}


a set of diffusion equations for the concentrations :math:`C_j`, for the substitutional elements
:math:`j = 2,\ldots, n-1`

.. math::
   :label: elphf:substitutional

   \underbrace{
	\frac{\partial C_j}{\partial t}
   }_{\text{transient}}
    &= \underbrace{
	D_{j}\nabla^2 C_j
	\vphantom{\frac{\partial C_j}{\partial t}}
    }_{\text{diffusion}} \\
    & \qquad + \underbrace{
	D_{j}\nabla\cdot 
	\frac{C_j}{1 - \sum_{\substack{k=2\\ k \neq j}}^{n-1} C_k}
	\left\{
	    \overbrace{
		\sum_{\substack{i=2\\ i \neq j}}^{n-1} \nabla C_i
	    }^{\text{counter diffusion}}
	   + 
	   \overbrace{ 
		C_n \left[
		    p'(\xi) \Delta\mu_{jn}^{\circ}
		    + g'(\xi) W_{jn}
		\right] \nabla\xi
		\vphantom{\sum_{\substack{i=2\\ i \neq j}}^{n-1} \nabla C_i}
	   }^{\text{phase transformation}}
	   +
	   \overbrace{
		C_n z_{jn} \nabla \phi
		\vphantom{\sum_{\substack{i=2\\ i \neq j}}^{n-1} \nabla C_i}
	   }^{\text{electromigration}}
	\right\}
    }_{\text{convection}}

a diffusion equation for the concentration :math:`C_{\text{e}^{-}}` of
electrons

.. math::
    :label: elphf:interstitial

    \underbrace{
	\frac{\partial C_{\text{e}^{-}}}{\partial t}
	\vphantom{\left\{
	    \overbrace{
		\left[p'(\xi)\right]
	    }^{\text{phase transformation}}
	\right\}}
    }_{\text{transient}}
    = \underbrace{
	D_{\text{e}^{-}}\nabla^2 C_{\text{e}^{-}} 
	\vphantom{\left\{
	    \overbrace{
		\left[p'(\xi)\right]
	    }^{\text{phase transformation}}
	\right\}}
    }_{\text{diffusion}} \\
    + \underbrace{
	D_{\text{e}^{-}}\nabla\cdot 
	C_{\text{e}^{-}}
	\left\{
	    \overbrace{
		\left[
		    p'(\xi) \Delta\mu_{\text{e}^{-}}^{\circ}
		    + g'(\xi) W_{\text{e}^{-}}
		\right] \nabla\xi
	    }^{\text{phase transformation}}
	    +
	    \overbrace{
		z_{\text{e}^{-}} \nabla \phi
		\vphantom{\left[p'(\xi)\right]}
	   }^{\text{electromigration}}
	\right\}
    }_{\text{convection}}

and Poisson's equation for the electrostatic potential :math:`\phi`

.. math::
    :label: elphf:Poisson

    \underbrace{
	\nabla\cdot\left(\epsilon\nabla\phi\right) 
	\vphantom{\sum_{j=1}^n z_j C_j}
    }_{\text{diffusion}}
    +
    \underbrace{
	\sum_{j=1}^n z_j C_j
    }_{\text{source}}
    = 0

:math:`M_\xi` is the phase field mobility, :math:`\kappa_\xi` is the phase
field gradient energy coefficient,
:math:`p'(\xi) = 30\xi^2\left(1-\xi\right)^2`, and
:math:`g'(\xi) = 2\xi\left(1-\xi\right)\left(1-2\xi\right)`.  For a given species
:math:`j`, :math:`\Delta\mu_j^{\circ}` is the standard chemical potential
difference between the electrode and electrolyte for a pure material,
:math:`W_j` is the magnitude of the energy barrier in the double-well
free energy function, :math:`z_j` is the valence, and :math:`D_{j}` is the
self diffusivity.  :math:`\Delta\mu_{jn}^{\circ}`, :math:`W_{jn}`, and
:math:`z_{jn}` are the differences of the respective quantities
:math:`\Delta\mu_{j}^{\circ}`, :math:`W_{j}`, and :math:`z_{j}` between
substitutional species :math:`j` and the solvent species :math:`n`.  The
total charge is denoted by :math:`\sum_{j=1}^n z_j C_j`.

Although unresolved stiffnesses make the full solution of this coupled
set of equations intractable in :term:`FiPy`, the following examples
demonstrate the setup and solution of various parts.
"""
