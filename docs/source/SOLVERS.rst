.. _SOLVERS:

=======
Solvers
=======

:term:`FiPy` requires either PETSc_, pyamgx_, Pysparse_, SciPy_, or
Trilinos_ solver suites to be installed in order to solve linear systems.
PETSc_ and Trilinos_ are the most complete of the
solvers due to their numerous preconditioning and solver capabilities and
they also allow :term:`FiPy` to :ref:`run in parallel <PARALLEL>`.
The :term:`Python` interface for PETSc_
is better maintained than for Trilinos_ and tends to be easier to install.
The sparse linear algebra solvers from the popular SciPy_ package are
widely available and easy to install. Although they do not perform as well
as the other suites and lack many of the features of PETSc_ or Trilinos_,
they may be the easiest linear solver choice when
first installing and testing :term:`FiPy`.
While the Pysparse_ linear solvers offer a modest advantage in serial, be
aware that they require :term:`Python` 2.7, which is no longer supported.
FiPy support for Pysparse_ will be dropped soon.
pyamgx_ offers the possibility
of solving sparse linear systems on the GPU; be aware that both
hardware and software configuration is non-trivial.

:term:`FiPy` chooses the solver suite based on system availability or based
on the user supplied :ref:`FlagsAndEnvironmentVariables`. For example,
passing ``--no-pysparse``::

    $ python -c "from fipy import *; print DefaultSolver" --no-pysparse
    <class 'fipy.solvers.trilinos.linearGMRESSolver.LinearGMRESSolver'>

uses a :ref:`TRILINOS` solver. Setting :envvar:`FIPY_SOLVERS`
to ``scipy``::

    $ FIPY_SOLVERS=scipy
    $ python -c "from fipy import *; print DefaultSolver"
    <class 'fipy.solvers.scipy.linearLUSolver.LinearLUSolver'>

uses a SciPy_ solver. Suite-specific solver classes can also
be imported and instantiated overriding any other directives. For
example::

    $ python -c "from fipy.solvers.scipy import DefaultSolver; \
    >   print DefaultSolver" --no-pysparse
    <class 'fipy.solvers.scipy.linearLUSolver.LinearLUSolver'>

uses a SciPy_ solver regardless of the command line
argument. In the absence of :ref:`FlagsAndEnvironmentVariables`,
:term:`FiPy`'s order of precedence when choosing the
solver suite for generic solvers is PySparse_ followed by
PETSc_, Trilinos_, SciPy_, PyAMG_, and pyamgx_.

.. note:: :term:`FiPy` has not been designed to mix different solver
   suites during a given problem run

Serial performance is compared for the different suites (see
:ref:`PARALLEL` for an analysis of parallel performance).

.. plot:: pyplots/scaling.py
   :align: center
   :alt: Wall time vs mesh size on a log-log plot.

   Comparison of serial performance for different solver suites, solvers
   and preconditioners [#Binary]_.  (a) Total elapsed time, (b) time to
   prepare the matrix, and (c) time to solve the matrix as a function of
   mesh size.  and different versions of :term:`FiPy` [#FIPYversion]_.  The
   light-colored bands in (a) and (b) denote plus-or-minus one standard
   deviation, based on five runs.  The fine traces in (c) denote each
   available combination of solver and preconditioner.  The thick traces in
   (c) show the ``LinearGMRESSolver`` and ``JacobiPreconditioner`` for each
   solver suite.

We can see:

- For sufficiently large problems, building the matrix can be expected to
  scale as the number of cells :math:`N` and solving the matrix should scale
  as :math:`N\,\ln N`.  There are not sufficient data points to
  differentiate these slopes.
- Below about 1000 cells, the time to prepare the matrix is insensitive to
  mesh size and this dominates the overall elapsed time.
- There is nearly three orders of magnitude between the fastest
  solver/preconditioner and the slowest.  This particular problem is not
  especially sensitive to choice of solver and preconditioner, as preparing
  the matrix takes 70 % to 80 % of the overall time, but it can be worth
  optimizing the choice for more complex systems of equations.
- `PETSc 3.19 <https://petsc.org/release/changes/319/>`_ introduced
  changes to "provide reasonable performance when no preallocation
  information is provided". Our experience is opposite that;
  :term:`FiPy` did not supply preallocation information prior to version
  4.0, but matrix preparation performance was fine with older :ref:`PETSC`
  releases.  Matrix preparation time is terrible when older :term:`FiPy` is
  combined with newer :ref:`PETSC`.  :term:`FiPy` 4.0 does supply preallocation
  information and matrix preparation time is comparable for all tested
  versions of :ref:`PETSC`.
- The :ref:`PETSC` and :ref:`PYSPARSE` solvers are about three times faster
  than the :ref:`SCIPY` or :ref:`TRILINOS` solvers (the :ref:`TRILINOS`
  solve time is generally good, but we have not been able to achieve the
  same level of performance for matrix preparation).

.. [#Binary] Calculations are of diffusion of a binary alloy in a frozen
   two-phase field.  Solutions are on a square
   :class:`~fipy.meshes.grid2D.Grid2D`.  All available solvers and
   preconditioners are attempted.  Solution tolerance is ``1e-10`` using
   the ``"RHS"`` :ref:`convergence criterion <CONVERGENCE>`.  Simulations
   were run on an AMD Epyc 7702 CPU with 64 cores featuring two-thread
   Simultaneous Multi-Threading (SMT) and 512 GB of memory.

.. [#FIPYversion] :term:`FiPy` version 3.4.4 has different interpretations
   of :ref:`CONVERGENCE` for different solver suites (and even for
   different solvers). Benchmarks used a patched version
   (`371d28468 <https://github.com/usnistgov/fipy/tree/371d28468>`_) that
   provided more logging information and normalized interpretation of
   tolerance, but without any of the improvements in matrix and solver
   efficiency of version 4.0.

.. _PETSC:

-----
PETSc
-----

https://www.mcs.anl.gov/petsc

:term:`PETSc` (the Portable, Extensible Toolkit for Scientific Computation)
is a suite of data structures and routines for the scalable (parallel)
solution of scientific applications modeled by partial differential
equations.  It employs the :term:`MPI` standard for all message-passing
communication (see :ref:`PARALLEL` for more details).

.. attention:: :term:`PETSc` requires the :term:`petsc4py` and :term:`mpi4py`
   interfaces.

.. note:: While, for consistency with other solver suites, :term:`FiPy` does
   implement some preconditioner objects for :term:`PETSc`, you can also
   simply pass one of the `PCType strings`_ in the `precon=` argument when
   declaring the solver.

.. _PCType strings: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html

.. _PYSPARSE:

--------
Pysparse
--------

http://pysparse.sourceforge.net

:term:`Pysparse` is a fast serial sparse matrix library for :term:`Python`.
It provides several sparse matrix storage formats and conversion methods.
It also implements a number of iterative solvers, preconditioners, and
interfaces to efficient factorization packages. The only requirement to
install and use :term:`Pysparse` is :term:`NumPy`.

.. warning::

   :term:`Pysparse` is archaic and limited to :ref:`RunningUnderPython2`.

.. warning::

   Support for :term:`Python` 2.7 and, thus, for :term:`Pysparse`
   will be dropped soon.

.. warning::

   :term:`FiPy` requires version 1.0 or higher of :term:`Pysparse`.

.. _SCIPY:

-----
SciPy
-----

http://www.scipy.org/

The :mod:`scipy.sparse` module provides a basic set of serial Krylov
solvers and a limited collection of preconditioners.

.. _PYAMG:

-----
PyAMG
-----

http://code.google.com/p/pyamg/

The :term:`PyAMG` package provides adaptive multigrid preconditioners that
can be used in conjunction with the SciPy_ solvers.

.. _PYAMGX:

------
pyamgx
------

https://pyamgx.readthedocs.io/

The :term:`pyamgx` package is a :term:`Python` interface to the 
NVIDIA `AMGX <https://github.com/NVIDIA/AMGX>`_
library.  :term:`pyamgx` can be used to construct complex solvers and
preconditioners to solve sparse sparse linear systems on the GPU.

.. _TRILINOS:

--------
Trilinos
--------

http://trilinos.sandia.gov

:term:`Trilinos` provides a more complete set of solvers and
preconditioners than either Pysparse_ or
SciPy_. :term:`Trilinos` preconditioning allows for iterative
solutions to some difficult problems that Pysparse_ and
SciPy_ cannot solve, and it enables parallel execution of
:term:`FiPy` (see :ref:`PARALLEL` for more details).

.. attention::

   Be sure to build or install the :term:`PyTrilinos` interface to
   :term:`Trilinos`.

.. attention::

   :term:`Trilinos` is a large software suite with its own set of
   prerequisites, and can be difficult to set up. It is not necessary
   for most problems, and is **not** recommended for a basic install
   of :term:`FiPy`.

.. attention::

   :term:`Trilinos` *must* be compiled with :term:`MPI` support for
   :ref:`PARALLEL`.

.. tip::

   :term:`Trilinos` parallel efficiency is somewhat improved by also
   installing :term:`Pysparse`.

.. note::

    :term:`Trilinos` solvers frequently give intermediate output that
    :term:`FiPy` cannot suppress. The most commonly encountered
    messages are

     ``Gen_Prolongator warning : Max eigen <= 0.0``
        which is not significant to :term:`FiPy`.

     ``Aztec status AZ_loss: loss of precision``
        which indicates that there was some difficulty in solving the
        problem to the requested tolerance due to precision limitations,
        but usually does not prevent the solver from finding an adequate
        solution.

     ``Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned``
        which indicates that GMRES is having trouble with the problem, and
        may indicate that trying a different solver or preconditioner may
        give more accurate results if GMRES fails.

     ``Aztec status AZ_breakdown: numerical breakdown``
        which usually indicates serious problems solving the equation which
        forced the solver to stop before reaching an adequate solution.
        Different solvers, different preconditioners, or a less restrictive
        tolerance may help.

.. _CONVERGENCE:

-----------
Convergence
-----------

Different solver suites take different approaches to testing convergence.
We endeavor to harmonize this behavior by allowing the strings in the
"criterion" column to be passed as an argument when instantiating a
:class:`~fipy.solvers.solver.Solver`.  Convergence is detected if
``residual < tolerance * scale``.

.. raw:: latex

    \begin{landscape}

.. csv-table:: Residual Criteria
   :file: _static/residuals.csv
   :widths: auto
   :header-rows: 1
   :stub-columns: 1
   :class: wideshow longtable

.. raw:: latex

    \end{landscape}

.. note:: PyAMG_ is a set of preconditioners applied on top of SciPy_, so
   is not explicitly included in these tables.



``default``
===========

The setting ``criterion="default"`` applies the same scaling (``RHS``) to
all solvers.  This behavior is new in version |release|; prior to that, the
default behavior was the same as ``criterion="legacy"``.

``legacy``
==========

The setting ``criterion="legacy"`` restores the behavior of FiPy prior to
version |release| and is equivalent to what the particular suite and solver
does if not specifically configured.  The ``legacy`` row of the table is a
best effort at documenting what will happen.

.. note::

    - All LU solvers use ``"initial"`` scaling.
    - PySparse_ has two different groups of solvers,
      with different scaling.
    - PETSc_ accepts |KSP_NORM_DEFAULT|_ in order to
      "use the default for the current ``KSPType``".  Discerning the actual
      behavior would require burning the code in a bowl of chicken entrails.
      (It is reasonable to assume |KSP_NORM_PRECONDITIONED|_ for
      left-preconditioned solvers and |KSP_NORM_UNPRECONDITIONED|_
      otherwise.
    - Even the PETSc_ documentation says that |KSP_NORM_NATURAL|_ is `"weird"
      <https://petsc.org/main/manualpages/KSP/KSPCGS/#developer-note>`_).

``absolute_tolerance``
======================

PETSc_ and SciPy_ Krylov solvers accept an additional
``absolute_tolerance`` parameter, such that convergence is detected if
``residual < max(tolerance * scale, absolute_tolerance``).

``divergence_tolerance``
========================

PETSc_ Krylov solvers accept a third ``divergence_tolerance`` parameter,
such that a divergence is detected if ``residual > divergence_tolerance *
scale``.  Because of `the way the convergence test is coded
<https://gitlab.com/petsc/petsc/-/blob/main/src/ksp/ksp/interface/iterativ.c#L1598>`_,
if the initial residual is much larger than the norm of the right-hand-side
vector, PETSc_ will abort with |KSP_DIVERGED_DTOL|_ without ever trying to
solve.  If this occurs, either ``divergence_tolerance`` should be increased
or another convergence criterion should be used.

.. note::

   See :mod:`examples.diffusion.mesh1D`,
   :mod:`examples.diffusion.steadyState.mesh1D.inputPeriodic`,
   :mod:`examples.elphf.diffusion.mesh1D`,
   :mod:`examples.elphf.phaseDiffusion`, :mod:`examples.phase.binary`,
   :mod:`examples.phase.quaternary`, and
   :mod:`examples.reactiveWetting.liquidVapor1D` for several examples where
   :code:`criterion="initial"` is used to address this situation.

.. note::

   ``divergence_tolerance`` never caused a problem in previous versions of
   :term:`FiPy` because the default behavior of PETSc_ is to zero out the
   initial guess before trying to solve and then never do a test against
   ``divergence_tolerance``.  This resulted in behavior (number of
   iterations and ultimate residual) that was very different from the other
   solver suites and so :term:`FiPy` now directs PETSc to use the initial
   guess.

Reporting
=========

Different solver suites also report different levels of detail about why
they succeed or fail.  This information is captured as a
:class:`~fipy.solvers.convergence.Convergence` or
:class:`~fipy.solvers.convergence.Divergence` property of the
:class:`~fipy.solvers.solver.Solver` after calling
:meth:`~fipy.terms.term.Term.solve` or
:meth:`~fipy.terms.term.Term.sweep`.

.. raw:: latex

    \begin{landscape}

.. tabularcolumns:: \Y{.25}\Y{.10}\Y{.22}\Y{.16}\Y{.09}\Y{.06}\Y{.12}

.. csv-table:: Convergence Status Codes
   :file: _static/solver_convergence.csv
   :widths: auto
   :header-rows: 1
   :stub-columns: 1
   :class: wideshow longtable

.. raw:: latex

    \end{landscape}


.. raw:: latex

    \begin{landscape}

.. tabularcolumns:: \Y{.25}\Y{.10}\Y{.22}\Y{.16}\Y{.09}\Y{.06}\Y{.12}

.. csv-table:: Divergence Status Codes
   :file: _static/solver_divergence.csv
   :widths: auto
   :header-rows: 1
   :stub-columns: 1
   :class: wideshow longtable

.. raw:: latex

    \end{landscape}

.. |KSP_NORM_UNPRECONDITIONED|  replace:: :literal:`KSP_NORM_UNPRECONDITIONED`
.. _KSP_NORM_UNPRECONDITIONED:  https://petsc.org/main/docs/manualpages/KSP/KSP_NORM_UNPRECONDITIONED/
.. |KSP_NORM_PRECONDITIONED|  replace:: :literal:`KSP_NORM_PRECONDITIONED`
.. _KSP_NORM_PRECONDITIONED:  https://petsc.org/main/docs/manualpages/KSP/KSP_NORM_PRECONDITIONED/
.. |KSP_NORM_NATURAL|  replace:: :literal:`KSP_NORM_NATURAL`
.. _KSP_NORM_NATURAL:  https://petsc.org/main/docs/manualpages/KSP/KSP_NORM_NATURAL/
.. |KSP_NORM_DEFAULT|  replace:: :literal:`KSP_NORM_DEFAULT`
.. _KSP_NORM_DEFAULT:  https://petsc.org/main/manualpages/KSP/KSPNormType/

.. [#KSP_Convergence_Tests] https://petsc.org/release/docs/manual/ksp/#sec-convergencetests

.. [#AMGX_convergence]   *AMGX REFERENCE MANUAL*: 2.3 General Settings: ``convergence``,
   October 2017, API Version 2,
   https://github.com/NVIDIA/AMGX/blob/main/doc/AMGX_Reference.pdf

.. [#SciPy_Convergence_Test]  https://github.com/scipy/scipy/blob/2d1d5b042a09e131ffe191726aa6829b33590970/scipy/sparse/linalg/_isolve/iterative.py#L30

.. [#AztecOO_convergence]  *AztecOO Users Guide*: 3.1  Aztec Options: ``options[AZ_conv]``,
   SAND REPORT SAND2004-3796, Updated August 2007,
   For AztecOO Version 3.6 in Trilinos Release 8.0,
   https://trilinos.github.io/pdfs/AztecOOUserGuide.pdf

.. [#FiPy_Convergence_Test] Implemented by :term:`FiPy` using intrinsic
   solver capabilities.

.. |KSP_CONVERGED_ITS|             replace:: :literal:`KSP_CONVERGED_ITS`
.. _KSP_CONVERGED_ITS:             https://petsc.org/main/docs/manualpages/KSP/KSP_CONVERGED_ITS/
.. |KSP_CONVERGED_ATOL|            replace:: :literal:`KSP_CONVERGED_ATOL`
.. _KSP_CONVERGED_ATOL:            https://petsc.org/main/docs/manualpages/KSP/KSP_CONVERGED_ATOL/
.. |KSP_CONVERGED_RTOL|            replace:: :literal:`KSP_CONVERGED_RTOL`
.. _KSP_CONVERGED_RTOL:            https://petsc.org/main/docs/manualpages/KSP/KSP_CONVERGED_RTOL/
.. |KSP_CONVERGED_ITERATING|       replace:: :literal:`KSP_CONVERGED_ITERATING`
.. _KSP_CONVERGED_ITERATING:       https://petsc.org/main/docs/manualpages/KSP/KSP_CONVERGED_ITERATING/
.. |KSP_DIVERGED_ITS|              replace:: :literal:`KSP_DIVERGED_ITS`
.. _KSP_DIVERGED_ITS:              https://petsc.org/main/docs/manualpages/KSP/KSP_DIVERGED_ITS/
.. |KSP_DIVERGED_PC_FAILED|        replace:: :literal:`KSP_DIVERGED_PC_FAILED`
.. _KSP_DIVERGED_PC_FAILED:        https://petsc.org/main/docs/manualpages/KSP/KSP_DIVERGED_PC_FAILED/
.. |KSP_DIVERGED_INDEFINITE_PC|    replace:: :literal:`KSP_DIVERGED_INDEFINITE_PC`
.. _KSP_DIVERGED_INDEFINITE_PC:    https://petsc.org/main/docs/manualpages/KSP/KSP_DIVERGED_INDEFINITE_PC/
.. |KSP_DIVERGED_INDEFINITE_MAT|   replace:: :literal:`KSP_DIVERGED_INDEFINITE_MAT`
.. _KSP_DIVERGED_INDEFINITE_MAT:   https://petsc.org/main/docs/manualpages/KSP/KSPConvergedReason/
.. |KSP_DIVERGED_NANORINF|         replace:: :literal:`KSP_DIVERGED_NANORINF`
.. _KSP_DIVERGED_NANORINF:         https://petsc.org/main/docs/manualpages/KSP/KSPConvergedReason/
.. |KSP_DIVERGED_BREAKDOWN|        replace:: :literal:`KSP_DIVERGED_BREAKDOWN`
.. _KSP_DIVERGED_BREAKDOWN:        https://petsc.org/main/docs/manualpages/KSP/KSP_DIVERGED_BREAKDOWN/
.. |KSP_DIVERGED_BREAKDOWN_BICG|   replace:: :literal:`KSP_DIVERGED_BREAKDOWN_BICG`
.. _KSP_DIVERGED_BREAKDOWN_BICG:   https://petsc.org/main/docs/manualpages/KSP/KSP_DIVERGED_BREAKDOWN_BICG/
.. |KSP_CONVERGED_HAPPY_BREAKDOWN| replace:: :literal:`KSP_CONVERGED_HAPPY_BREAKDOWN`
.. _KSP_CONVERGED_HAPPY_BREAKDOWN: https://petsc.org/main/docs/manualpages/KSP/KSPConvergedReason/
.. |KSP_DIVERGED_NULL|             replace:: :literal:`KSP_DIVERGED_NULL`
.. _KSP_DIVERGED_NULL:             https://petsc.org/main/docs/manualpages/KSP/KSPConvergedReason/
.. |KSP_DIVERGED_DTOL|             replace:: :literal:`KSP_DIVERGED_DTOL`
.. _KSP_DIVERGED_DTOL:             https://petsc.org/main/docs/manualpages/KSP/KSP_DIVERGED_DTOL/
.. |KSP_DIVERGED_NONSYMMETRIC|     replace:: :literal:`KSP_DIVERGED_NONSYMMETRIC`
.. _KSP_DIVERGED_NONSYMMETRIC:     https://petsc.org/main/docs/manualpages/KSP/KSP_DIVERGED_NONSYMMETRIC/

.. |AMGX_SOLVE_SUCCESS|            replace:: :literal:`AMGX_SOLVE_SUCCESS`
.. _AMGX_SOLVE_SUCCESS:            https://github.com/NVIDIA/AMGX/blob/main/doc/AMGX_Reference.pdf
.. |AMGX_SOLVE_FAILED|             replace:: :literal:`AMGX_SOLVE_FAILED`
.. _AMGX_SOLVE_FAILED:             https://github.com/NVIDIA/AMGX/blob/main/doc/AMGX_Reference.pdf
.. |AMGX_SOLVE_DIVERGED|           replace:: :literal:`AMGX_SOLVE_DIVERGED`
.. _AMGX_SOLVE_DIVERGED:           https://github.com/NVIDIA/AMGX/blob/main/doc/AMGX_Reference.pdf

.. |PySparse_2|                    replace:: :literal:`2`
.. _PySparse_2:                    http://pysparse.sourceforge.net/itsolvers.html
.. |PySparse_1|                    replace:: :literal:`1`
.. _PySparse_1:                    http://pysparse.sourceforge.net/itsolvers.html
.. |PySparse_0|                    replace:: :literal:`0`
.. _PySparse_0:                    http://pysparse.sourceforge.net/itsolvers.html
.. |PySparse_neg1|                 replace:: :literal:`-1`
.. _PySparse_neg1:                 http://pysparse.sourceforge.net/itsolvers.html
.. |PySparse_neg2|                 replace:: :literal:`-2`
.. _PySparse_neg2:                 http://pysparse.sourceforge.net/itsolvers.html
.. |PySparse_neg3|                 replace:: :literal:`-3`
.. _PySparse_neg3:                 http://pysparse.sourceforge.net/itsolvers.html
.. |PySparse_neg4|                 replace:: :literal:`-4`
.. _PySparse_neg4:                 http://pysparse.sourceforge.net/itsolvers.html
.. |PySparse_neg5|                 replace:: :literal:`-5`
.. _PySparse_neg5:                 http://pysparse.sourceforge.net/itsolvers.html
.. |PySparse_neg6|                 replace:: :literal:`-6`
.. _PySparse_neg6:                 http://pysparse.sourceforge.net/itsolvers.html

.. |SciPy_0|                       replace:: :literal:`0`
.. _SciPy_0:                       https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html
.. |SciPy_lt0|                     replace:: :literal:`<0`
.. _SciPy_lt0:                     https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html
.. |SciPy_gt0|                     replace:: :literal:`>0`
.. _SciPy_gt0:                     https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html

.. |AZ_normal|                     replace:: :literal:`AZ_normal`
.. _AZ_normal:                     https://trilinos.github.io/pdfs/AztecOOUserGuide.pdf
.. |AZ_maxits|                     replace:: :literal:`AZ_maxits`
.. _AZ_maxits:                     https://trilinos.github.io/pdfs/AztecOOUserGuide.pdf
.. |AZ_ill_cond|                   replace:: :literal:`AZ_ill_cond`
.. _AZ_ill_cond:                   https://trilinos.github.io/pdfs/AztecOOUserGuide.pdf
.. |AZ_breakdown|                  replace:: :literal:`AZ_breakdown`
.. _AZ_breakdown:                  https://trilinos.github.io/pdfs/AztecOOUserGuide.pdf
.. |AZ_loss|                       replace:: :literal:`AZ_loss`
.. _AZ_loss:                       https://trilinos.github.io/pdfs/AztecOOUserGuide.pdf
