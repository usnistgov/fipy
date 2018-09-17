.. _sec:NumericalSchemes:

Numerical Schemes
-----------------

The coefficients of equation :eq:`num:dap` must remain positive,
since an increase in a neighboring value must result in an increase in
:math:`\phi_P` to obtain physically realistic solutions.  Thus, the
inequalities :math:`a_A > 0` and :math:`a_A + F_f>0` must be satisfied.  The
|Peclet| number :math:`P_f \equiv F_f / D_f` is the ratio between convective
strength and diffusive conductance.  To achieve physically realistic
solutions, the inequality

.. math::
  :label: eqn:num:inq

   \frac{1}{1-\alpha_f} > P_f > -\frac{1}{\alpha_f}

must be satisfied.
The parameter :math:`\alpha_f` is defined by the chosen scheme, depending
on Equation :eq:`eqn:num:inq`. The various
differencing schemes are:

the central differencing scheme,
  where

  .. math::
     :label: eqn:num:cds

     \alpha_f = \frac{1}{2}

  so that :math:`|P_f|<2` satisfies Equation :eq:`eqn:num:inq`. Thus, the
  central differencing scheme is only numerically stable for a low
  values of :math:`P_f`.

the upwind scheme,
  where

  .. math::
     :label: eqn:num:ups

     \alpha_f = \begin{cases}
     1 & \text{if $P_f > 0$,} \\
     0 & \text{if $P_f < 0$.}
     \end{cases}

  Equation :eq:`eqn:num:ups` satisfies the inequality in
  Equation :eq:`eqn:num:inq` for all values of :math:`P_f`.  However the
  solution over predicts the diffusive term leading to excessive
  numerical smearing ("false diffusion").

the exponential scheme,
  where

  .. math::
     :label: eqn:num:exs

     \alpha_f = \frac{(P_f-1)\exp{(P_f)}+1}{P_f(\exp{(P_f)}-1)}.

  This formulation can be derived from the exact solution, and thus,
  guarantees positive coefficients while not over-predicting the
  diffusive terms. However, the computation of exponentials is slow and
  therefore a faster scheme is generally used, especially in higher
  dimensions.

the hybrid scheme,
  where

  .. math::
     :label: eqn:num:hys

     \alpha_f =
     \begin{cases}
         \frac{P_f-1}{P_f} & \text{if $P_f > 2$,} \\
         \frac{1}{2} & \text{if $|P_f| < 2$,} \\
         -\frac{1}{P_f} & \text{if $P_f < -2$.}
     \end{cases}

  The hybrid scheme is formulated by allowing :math:`P_f \rightarrow \infty`,
  :math:`P_f \rightarrow 0` and :math:`P_f \rightarrow -\infty` in the exponential
  scheme.  The hybrid scheme is an improvement on the upwind scheme,
  however, it deviates from the exponential scheme at :math:`|P_f|=2`.

the power law scheme,
  where

  .. math::
     :label: eqn:num:pls

     \alpha_f =
     \begin{cases}
         \frac{P_f-1}{P_f} & \text{if $P_f > 10$,} \\
         \frac{(P_f-1)+(1-P_f/10)^5}{P_f} & \text{if $0 < P_f < 10$,} \\
         \frac{(1-P_f/10)^5 - 1}{P_f} & \text{if $-10 < P_f < 0$,} \\
         -\frac{1}{P_f} & \text{if $P_f < -10$.}
     \end{cases}

  The power law scheme overcomes the inaccuracies of the hybrid scheme,
  while improving on the computational time for the exponential scheme.

.. warning::

   :class:`~fipy.terms.vanLeerConvectionTerm.VanLeerConvectionTerm` not
   mentioned and no discussion of explicit forms.

All of the numerical schemes presented here are available in :term:`FiPy`
and can be selected by the user.

.. |Peclet| unicode:: P U+00E9 clet
