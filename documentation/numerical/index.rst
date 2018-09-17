.. _chap:Numerics:

====================================
Theoretical and Numerical Background
====================================

This chapter describes the numerical methods used to solve equations
in the :term:`FiPy` programming environment. :term:`FiPy` uses the finite volume
method (FVM) to solve coupled sets of partial differential equations
(PDEs). For a good introduction to the FVM see Nick Croft's PhD
thesis :cite:`croftphd`, Patankar :cite:`patankar` or Versteeg and
Malalasekera :cite:`versteegMalalasekera`.

Essentially, the FVM consists of dividing the solution domain into
discrete finite volumes over which the state variables are
approximated with linear or higher order interpolations. The
derivatives in each term of the equation are satisfied with simple
approximate interpolations in a process known as discretization. The
(FVM) is a popular discretization technique employed to solve coupled
PDEs used in many application areas (*e.g.*, Fluid Dynamics).

The FVM can be thought of as a subset of the Finite Element Method
(FEM), just as the Finite Difference Method (FDM) is a subset of the FVM. A
system of equations fully equivalent to the FVM can be obtained with
the FEM using as weighting functions the characteristic functions of
FV cells, i.e., functions equal to unity :cite:`Mattiussi:1997`. Analogously,
the discretization of equations with the FVM reduces to the FDM on
Cartesian grids.

.. toctree::

   equation
   discret
   scheme
