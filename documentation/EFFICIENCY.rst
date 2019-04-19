.. _chap:Efficiency:

==========
Efficiency
==========

This section will present results and discussion of efficiency
evaluations with :term:`FiPy`. Programming in :term:`Python` allows greater
efficiency when designing and implementing new code, but it has some
intrinsic inefficiencies during execution as compared with the C or
FORTRAN programming languages. These inefficiencies can be minimized
by translating sections of code that are used frequently into C.

:term:`FiPy` has been tested against an in-house phase field code, written
at NIST, to model grain growth and subsequent impingement. This
problem can be executed by running::

    $ examples/phase/impingement/mesh20x20.py \
    > --numberOfElements=10000 --numberOfSteps=1000

from the base :term:`FiPy` directory. The in-house code was written by Ryo
Kobayashi and is used to generate the results presented in
:cite:`WarrenPolycrystal`.

The raw CPU execution times for 10 time steps are presented in the
following table. The run times are in seconds and the memory usage is
in kilobytes. The Kobayashi code is given the heading of FORTRAN while
:term:`FiPy` is run with and without inlining. The memory usage is for
:term:`FiPy` simulations with the :option:`--inline`. The :option:`--no-cache` flag is
on in all cases for the following table.

========== ============ ======================= ============= ============= =============
 Elements   FiPy (s)     FiPy                    FORTRAN (s)   FiPy          FORTRAN
                         :option:`--inline` (s)                memory (KiB)  memory (KiB)
========== ============ ======================= ============= ============= =============
    100       0.77        0.30                   0.0009         39316          772
    400       0.87        0.37                   0.0031         39664          828
   1600       1.4         0.65                   0.017          40656         1044
   6400       3.7         2.0                    0.19           46124         1880
  25600      19          10                      1.3            60840         5188
 102400      79          43                      4.6           145820        18436
========== ============ ======================= ============= ============= =============

The plain :term:`Python` version of :term:`FiPy`, which uses ``Numeric`` for all
array operations, is around 17 times slower than the FORTRAN
code. Using the :option:`--inline` flag, this penalty is reduced to about 9
times slower.

It is hoped that in future releases of :term:`FiPy` the process of C
inlining for ``Variable`` objects will be automated. This may result
in some efficiency gains, greater than we are seeing for this
particular problem since all the ``Variable`` objects will be
inlined. Recent analysis has shown that a ``Variable`` with multiple
operations could be up to 6 times faster at calculating its value when
inlined.

As presented in the above table, memory usage was also recorded for
each :term:`FiPy` simulation. From the table, once base memory usage is
subtracted, each cell requires approximately 1.4 kilobytes of memory.
The measurement of the maximum memory spike is hard with dynamic
memory allocation, so these figures should only be used as a very
rough guide. The FORTRAN memory usage is exact since memory is not
allocated dynamically.


Efficiency comparison between :option:`--no-cache` and :option:`--cache` flags
==============================================================================

This table shows results for efficiency tests when using the caching
flags. Examples with more variables involved in complex expressions
show the largest improvement in memory usage. The :option:`--no-cache`
option mainly prevents intermediate variables created due to binary
operations from caching their values. This results in large memory
gains while not effecting run times substantially. The table below is
with :option:`--inline` switched on and with 102400 elements for each case.
The :option:`--no-cache` flag is the default option.

.. currentmodule:: examples

======================================================== ==================== ================= ====================== ==================
 Example                                                  time per step        time per step      memory per cell       memory per cell
                                                          ``--no-cache`` (s)   ``--cache`` (s)    ``--no-cache`` (KiB)  ``--cache`` (KiB)
======================================================== ==================== ================= ====================== ==================
 :mod:`examples.phase.impingement.mesh20x20`                   4.3                  4.1               1.4                   2.3
 :mod:`examples.phase.anisotropy`                              3.5                  3.2               1.1                   1.9
 :mod:`examples.cahnHilliard.mesh2D`                           3.0                  2.5               1.1                   1.4
 :mod:`examples.levelSet.electroChem.simpleTrenchSystem`      62                   62                 2.0                   2.8
======================================================== ==================== ================= ====================== ==================

Efficiency discussion of Pysparse and Trilinos
==================================================================

Trilinos provides multigrid capabilities which are beneficial for
some problems, but has significant overhead compared to Pysparse.
The matrix-building step takes significantly longer in Trilinos,
and the solvers also have more overhead costs in memory and
performance than the equivalent Pysparse solvers. However, the
multigrid preconditioning capabilities of Trilinos can, in some
cases, provide enough of a speedup in the solution step to make up
for the overhead costs. This depends greatly on the specifics of
the problem, but is most likely in the cases when the problem is
large and when Pysparse cannot solve the problem with an iterative
solver and must use an LU solver, while Trilinos can still have
success with an iterative method.
