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
PETSc_, Trilinos_, SciPy_, and pyamgx_.

.. .. raw:: latex

    \begin{landscape}

.. .. tabularcolumns:: \Y{.25}\Y{.10}\Y{.22}\Y{.16}\Y{.09}\Y{.06}\Y{.12}

.. csv-table:: Solver Suite Features
   :file: _static/solver_features.csv
   :widths: 22, 13, 13, 13, 13, 13, 13
   :header-rows: 1
   :stub-columns: 1
   :class: wideshow longtable

..    :widths: 13,11,24,12,8,11,10,12

.. .. raw:: latex

    \end{landscape}

.. _mpi4py:             https://mpi4py.readthedocs.io/
.. |mpi4py|             replace:: :literal:`mpi4py`
.. _mpi4py:             https://mpi4py.readthedocs.io/
.. |petsc4py|           replace:: :literal:`petsc4py`
.. _petsc4py:           https://petsc4py.readthedocs.io/
.. |pyamg.|             replace:: :literal:`pyamg`
.. _pyamg.:             https://pyamg.readthedocs.io/
.. |pyamgx.|            replace:: :literal:`pyamgx`
.. _pyamgx.:            https://pyamgx.readthedocs.io/
.. |pysparse.|          replace:: :literal:`pysparse`
.. _pysparse.:          https://pysparse.sourceforge.net/
.. |PyTrilinos|         replace:: :literal:`PyTrilinos`
.. _PyTrilinos:         https://trilinos.github.io/pytrilinos.html
.. |PyTrilinos2.|       replace:: :literal:`PyTrilinos2`
.. _PyTrilinos2.:       https://trilinos.github.io/pytrilinos2.html
.. |scipy.|             replace:: :literal:`scipy`
.. _scipy.:             https://docs.scipy.org/doc/scipy/

.. [#pyamgx_parallel] While AMGX_ matrix solve takes advantage of GPU
   parallelism, the ``pyamgx`` library uses :ref:`SCIPY` to build the matrix
   and thus suffers a significant serial bottleneck.

.. [#37_38]  |PyTrilinos|_ may be compatible with newer versions of Python,
   but these are the most recent versions we've been able to get to install
   using :term:`conda` (3.7 on :term:`linux` and 3.8 on :term:`macOS`).

.. [#PyTrilinos2]  There is a more actively developed |PyTrilinos2.|_
   package, which may be compable with more recent versions of
   :term:`Python`, but :term:`FiPy` does not yet work with it.

.. [#PySparse4Trilinos] Trilinos parallel efficiency is somewhat improved
   by also installing |pysparse.|_.

.. note:: :term:`FiPy` has not been designed to mix different solver
   suites during a given problem run

.. _PETSC:

-----
PETSc
-----

https://www.mcs.anl.gov/petsc

PETSc (the Portable, Extensible Toolkit for Scientific Computation)
is a suite of data structures and routines for the scalable (parallel)
solution of scientific applications modeled by partial differential
equations.  It employs the :term:`MPI` standard for all message-passing
communication (see :ref:`PARALLEL` for more details).

.. _PYAMG:

-----
PyAMG
-----

https://pyamg.readthedocs.io/

The PyAMG package provides adaptive multigrid preconditioners that
can be used in conjunction with the SciPy_ solvers. While PyAMG also has
solvers, they are not currently implemented in :term:`FiPy`.

.. _PYAMGX:

------
pyamgx
------

https://pyamgx.readthedocs.io/

The pyamgx package is a :term:`Python` interface to the NVIDIA AMGX_
library.  pyamgx can be used to construct complex solvers and
preconditioners to solve sparse sparse linear systems on the GPU.

.. _AMGX: https://github.com/NVIDIA/AMGX

.. _PYSPARSE:

--------
Pysparse
--------

http://pysparse.sourceforge.net

Pysparse is a fast serial sparse matrix library for :term:`Python`.
It provides several sparse matrix storage formats and conversion methods.
It also implements a number of iterative solvers, preconditioners, and
interfaces to efficient factorization packages. The only requirement to
install and use Pysparse is :term:`NumPy`.

.. warning::

   Pysparse is archaic and limited to :ref:`RunningUnderPython2`.
   Support for :term:`Python` 2.7 and, thus, for Pysparse
   will be dropped soon.

.. _SCIPY:

-----
SciPy
-----

https://docs.scipy.org/doc/scipy/reference/sparse.html

The :mod:`scipy.sparse` module provides a basic set of serial Krylov
solvers and a limited collection of preconditioners.

.. _TRILINOS:

--------
Trilinos
--------

https://trilinos.github.io/

Trilinos provides a complete set of sparse matrices, solvers, and
preconditioners.  Trilinos preconditioning allows for iterative solutions
to some difficult problems, and it enables parallel execution of
:term:`FiPy` (see :ref:`PARALLEL` for more details).

----------------------
Performance Comparison
----------------------

Comparing different solver suites, or even different solvers, has
historically been difficult.  The different suites have different
interpretations of :ref:`CONVERGENCE` and tolerance.  :term:`FiPy` 4.0
harmonizes the different suites so that, to the greatest extent possible,
all interpret :ref:`CONVERGENCE` and tolerance the same way.  In the course
of doing this, a number of inefficiencies were found in the way that
:term:`FiPy` built sparse matrices.  To see the impact of these changes, we
examine the serial and parallel scaling performance of the different solver
suites for two different benchmark problems.

Serial Performance
==================

Serial performance is compared for the different suites.

The following plot shows the serial scaling behavior for the different
solvers.  We compare solution time vs number mesh cells for a diffusion
problem.

.. plot:: pyplots/serial_scaling.py
   :align: center
   :alt: Wall time vs mesh size on a log-log plot.

   Comparison of serial performance for different solver suites, solvers
   and preconditioners, and different versions of :term:`FiPy`
   [#FIPYversion]_.  (a) Total elapsed time, (b) time to prepare the
   matrix, and (c) time to solve the matrix as functions of mesh size.
   [#Binary]_

We can see:

- For sufficiently large problems, building the matrix can be expected to
  scale as the number of cells :math:`N` and solving the matrix should scale
  as :math:`N\,\ln N`.  There are not enough data points to
  differentiate these slopes.
- Below about 1000 cells, the time to prepare the matrix is insensitive to
  mesh size and this dominates the overall elapsed time.
- There is nearly three orders of magnitude between the fastest
  solver/preconditioner and the slowest.  This particular problem is not
  especially sensitive to choice of solver and preconditioner, as preparing
  the matrix takes the majority of the overall time, but it can be worth
  optimizing the choice for more complex systems of equations.
- Matrix preparation time is terrible when older :term:`FiPy` is
  combined with newer :ref:`PETSC`.  `PETSc 3.19
  <https://petsc.org/release/changes/319/>`_ introduced changes to "provide
  reasonable performance when no preallocation information is provided".
  Our experience is opposite that; :term:`FiPy` did not supply
  preallocation information prior to version 4.0, but matrix preparation
  performance was fine with older :ref:`PETSC` releases.  :term:`FiPy` 4.0
  does supply preallocation information and matrix preparation time is
  comparable for all tested versions of :ref:`PETSC`.
- There is considerable dispersion about the mean solve time for different
  solvers and preconditioners.  On the other hand, the time to prepare the
  matrix is insensitive to the choice of solver and preconditioner and
  shows a high degree of consistency from run to run.

.. plot:: pyplots/serial_fraction.py
   :align: center
   :alt: Fraction of time spent preparing matrix vs mesh size on a linear-log plot

   Ratio of time to prepare the matrix to the combined time to prepare and
   solve the matrix for different solver suites, solvers
   and preconditioners, and different versions of :term:`FiPy`
   [#FIPYversion]_ [#Binary]_.  The thick lines highlight ``LinearCGSolver``
   with no preconditioner, one of the better-performing combinations
   available in all suites.

In principle, we'd like to spend as little time preparing the matrix,
relative to solving it, as possible.  This metric can be deceptive.  If we
compare the case of unpreconditioned ``LinearCGSolver``, one of the fastest
combinations for all suites *for this problem*, we see that :ref:`Trilinos`
has the lowest ratio of prepare to elapsed time.  Examination of elapsed
time, the quantity we really care about, shows that :ref:`Trilinos` takes
three times as long to both prepare and solve as :ref:`PySparse` or
:ref:`SciPy` and twice as long as :ref:`PETSc`.

For your own work, focus on identifying the solver and preconditioner with
the lowest overall time to build and solve; this will counterintuitively
have the highest ratio of prepare-to-elapsed time.  Prepare time to elapsed
time is a more useful metric for the :term:`FiPy` developers; just as
:term:`FiPy` 4.0 brought considerable reductions in matrix build time, we
will continue to seek opportunities to optimize.

Parallel Performance
====================

The following plot shows the scaling behavior for multiple
processors.  We compare solution time vs number of Slurm_ tasks (available
cores) for a `Method of Manufactured Solutions Allen-Cahn problem`_.

.. plot:: pyplots/parallel_scaling.py
   :align: center
   :alt: "Speedup" relative to PySparse versus number of tasks (processes) on a log-log plot.

   Parallel scaling behavior of different solver packages and different
   versions of :term:`FiPy` [#FIPYversion]_ [#MMS]_.

A few things can be observed in this plot:

- :ref:`PETSc`, :ref:`PySparse`, :ref:`Trilinos`, and :ref:`SciPy` have
  comparable serial performance, with :ref:`SciPy` edging out the other
  three for this particular problem.

- :term:`FiPy` 4.0 is roughly the same speed in serial, but more than
  twice as fast in parallel compared to :term:`FiPy` 3.4.4 when using the
  :ref:`PETSC` solvers.  :term:`FiPy` 4.0 is roughly twice as fast
  using the :ref:`TRILINOS` solvers, whether in serial or parallel.

- :term:`FiPy` 4.0
  exhibits better parallel scaling than :term:`FiPy` 3.4.4.  `Amdahl's
  Law`_, :math:`\text{speedup} = p / (1 + \sigma(p - 1))`, does not fit the
  performance data nearly as well as `Gunther's Law`_,
  :math:`\text{speedup} = p / (1 + \sigma(p - 1) + \kappa p (p-1))`, where
  :math:`p` is the number of parallel tasks, :math:`\sigma` is the fraction
  limited by serial processes, and :math:`\kappa` is `"coherency" (which is
  somewhat nebulous)`_.

  .. table:: Parallel scaling fitting parameters (smaller numbers are better)

     +------------+----------+------------+------------+------------+
     |            |          | Amdahl     |         Gunther         |
     +------------+----------+------------+------------+------------+
     |            |          | serial / % | serial / % | coherency  |
     +============+==========+============+============+============+
     | FiPy 3.4.4 | petsc    | 4.7(3)     | 0.91(9)    | 0.00078(2) |
     +            +----------+------------+------------+------------+
     |            | trilinos | 2.6(1)     | 0.8(1)     | 0.00034(2) |
     +------------+----------+------------+------------+------------+
     | FiPy 4.0   | petsc    | 1.70(8)    | 0.13(7)    | 0.00028(1) |
     +            +----------+------------+------------+------------+
     |            | trilinos | 2.2(1)     | 0.4(1)     | 0.00032(3) |
     +------------+----------+------------+------------+------------+


At least one source of less-than-optimal scaling is that our
"``...Grid...``" meshes parallelize by dividing the mesh into slabs, which
leads to more communication overhead than more compact partitions.  The
"``...Gmsh...``" meshes partition more efficiently, but carry more overhead
in other ways.  We'll be making efforts to improve the partitioning of the
"``...Grid...``" meshes in a future release.

These results are likely both problem and architecture dependent.  You
should develop an understanding of the scaling behavior of your own codes
before doing "production" runs.

.. _Method of Manufactured Solutions Allen-Cahn problem:  https://pages.nist.gov/pfhub/benchmarks/benchmark7.ipynb
.. _Slurm: https://slurm.schedmd.com
.. _Windows Subsystem for Linux: https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux
.. _Amdahl's Law: https://en.wikipedia.org/wiki/Amdahl%27s_law
.. _Gunther's Law: https://doi.org/10.48550/arXiv.0808.1431
.. _"coherency" (which is somewhat nebulous): https://learn.microsoft.com/en-us/archive/blogs/ddperf/parallel-scalability-isnt-childs-play-part-2-amdahls-law-vs-gunthers-law

.. [#FIPYversion] :term:`FiPy` version 3.4.4 has different interpretations
   of :ref:`CONVERGENCE` for different solver suites (and even for
   different solvers). Benchmarks used a patched version
   (`371d28468 <https://github.com/usnistgov/fipy/tree/371d28468>`_) that
   provided more logging information and normalized interpretation of
   tolerance, but without any of the improvements in matrix and solver
   efficiency of version 4.0.

.. [#Binary] Calculations are of diffusion of a binary alloy in a frozen
   two-phase field.  Solutions are on a square
   :class:`~fipy.meshes.grid2D.Grid2D`.  The initial condition is sampled
   from the center of a well-evolved :math:`1024\times 1024`
   `nucleation simulation
   <https://pages.nist.gov/pfhub/benchmarks/benchmark8.ipynb/>`_.
   All available solvers and
   preconditioners are attempted.  Solution tolerance is ``1e-10`` using
   the ``"RHS"`` :ref:`convergence criterion <CONVERGENCE>`.  Simulations
   were run on an AMD Epyc 7702 CPU with 64 cores featuring two-thread
   Simultaneous Multi-Threading (SMT) and 512 GB of memory.
   :ref:`OMP_NUM_THREADS was set to 1 <THREADS_VS_RANKS>`.

.. [#MMS] Calculations are of a
   `Method of Manufactured Solutions Allen-Cahn problem`_.  Solutions are
   on a :math:`2048\times 1024` :class:`~fipy.meshes.grid2D.Grid2D`
   and the ``LinearCGSolver`` with no preconditioner is used for
   all solver suites.  Solution tolerance is ``1e-10`` using the ``"RHS"``
   :ref:`convergence criterion <CONVERGENCE>`.  Five replicates of each
   simulation were run on an AMD Epyc 7702 CPU with 64 cores featuring
   two-thread Simultaneous Multi-Threading (SMT) and 512 GB of memory.
   :ref:`OMP_NUM_THREADS was set to 1 <THREADS_VS_RANKS>`.

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
all solvers.  This behavior is new in :term:`FiPy` 4.0; prior to that, the
default behavior was the same as ``criterion="legacy"``.

``legacy``
==========

The setting ``criterion="legacy"`` restores the behavior of :term:`FiPy`
prior to version 4.0 and is equivalent to what the particular suite and solver
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
      otherwise.)
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
