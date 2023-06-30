.. _SOLVERS:

=======
Solvers
=======

:term:`FiPy` requires either :term:`Pysparse`, :term:`SciPy` or
:term:`Trilinos` to be installed in order to solve linear systems.
From our experiences, :term:`FiPy` runs most efficiently in serial
when :term:`Pysparse` is the linear solver. :term:`Trilinos` is the
most complete of the three solvers due to its numerous preconditioning
and solver capabilities and it also allows :term:`FiPy` to :ref:`run
in parallel <PARALLEL>`. Although less efficient than :term:`Pysparse`
and less capable than :term:`Trilinos`, :term:`SciPy` is a very
popular package, widely available and easy to install. For this
reason, :term:`SciPy` may be the best linear solver choice when first
installing and testing :term:`FiPy` (and it is the only viable solver
under `Python 3.x`_).

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

uses a :ref:`SCIPY` solver. Suite-specific solver classes can also
be imported and instantiated overriding any other directives. For
example::

    $ python -c "from fipy.solvers.scipy import DefaultSolver; \
    >   print DefaultSolver" --no-pysparse
    <class 'fipy.solvers.scipy.linearLUSolver.LinearLUSolver'>

uses a :ref:`SCIPY` solver regardless of the command line
argument. In the absence of :ref:`FlagsAndEnvironmentVariables`,
:term:`FiPy`'s order of precedence when choosing the
solver suite for generic solvers is :ref:`PYSPARSE` followed by
:ref:`TRILINOS`, :ref:`PYAMG` and :ref:`SCIPY`.

.. _Python 3.x:   http://docs.python.org/py3k/

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

.. note:: :term:`FiPy` does not implement any precoditioner objects for
   :term:`PETSc`. Simply pass one of the `PCType strings`_ in the
   `precon=` argument when declaring the solver.

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
install and use Pysparse is :term:`NumPy`.

.. warning::

   :term:`FiPy` requires version 1.0 or higher of :term:`Pysparse`.

.. _SCIPY:

-----
SciPy
-----

http://www.scipy.org/

The :mod:`scipy.sparse` module provides a basic set of serial Krylov
solvers, but no preconditioners.

.. _PYAMG:

-----
PyAMG
-----

http://code.google.com/p/pyamg/

The :term:`PyAMG` package provides adaptive multigrid preconditioners that
can be used in conjunction with the :term:`SciPy` solvers.

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
preconditioners than either :term:`Pysparse` or
:term:`SciPy`. :term:`Trilinos` preconditioning allows for iterative
solutions to some difficult problems that :term:`Pysparse` and
:term:`SciPy` cannot solve, and it enables parallel execution of
:term:`FiPy` (see :ref:`PARALLEL` for more details).

.. attention::

   Be sure to build or install the :term:`PyTrilinos` interface to
   :term:`Trilinos`.

.. attention::

   :term:`FiPy` runs more efficiently when :term:`Pysparse` is
   installed alongside :term:`Trilinos`.

.. attention::

   :term:`Trilinos` is a large software suite with its own set of
   prerequisites, and can be difficult to set up. It is not necessary
   for most problems, and is **not** recommended for a basic install
   of :term:`FiPy`.

.. attention::

   :term:`Trilinos` *must* be compiled with :term:`MPI` support for
   :ref:`PARALLEL`.

.. tip::

   :term:`Trilinos` parallel efficiency is greatly improved by also
   installing :term:`Pysparse`.  If :term:`Pysparse` is not installed, be
   sure to use the ``--no-pysparse`` flag.

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
