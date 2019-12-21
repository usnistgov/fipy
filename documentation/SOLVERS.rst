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
equations.  It employs the MPI standard for all message-passing
communication (see :ref:`PARALLEL` for more details).

.. attention:: :term:`PETSc` requires the :term:`petsc4py` and :ref:`mpi4py`
   interfaces.

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

:term:`Trilinos` requires `cmake <http://www.cmake.org/>`_, :term:`NumPy`,
and `swig <http://www.swig.org/>`_. The following are the minimal steps to
build and install :term:`Trilinos` (with :term:`PyTrilinos`) for
:term:`FiPy`::

    $ cd trilinos-X.Y/
    $ SOURCE_DIR=`pwd`
    $ mkdir BUILD_DIR
    $ cd BUILD_DIR
    $ cmake \
    >   -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    >   -D Trilinos_ENABLE_PyTrilinos:BOOL=ON \
    >   -D BUILD_SHARED_LIBS:BOOL=ON \
    >   -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
    >   -D TPL_ENABLE_MPI:BOOL=ON \
    >   -D Trilinos_ENABLE_TESTS:BOOL=ON \
    >   -D DART_TESTING_TIMEOUT:STRING=600 \
    >   ${SOURCE_DIR}
    $ make
    $ make install

Depending on your platform, other options may be helpful or necessary;
see the :term:`Trilinos` user guide available from
http://trilinos.sandia.gov/documentation.html, or
http://trilinos.sandia.gov/packages/pytrilinos/faq.html for more
in-depth documentation.

.. note::

    Trilinos can be installed in a non-standard location by adding the
    :samp:`-D CMAKE_INSTALL_PREFIX:PATH=${{INSTALL_DIR}}` and
    :samp:`-D PyTrilinos_INSTALL_PREFIX:PATH=${{INSTALL_DIR}}` flags
    to the configure step. If :term:`Trilinos` is installed in a
    non-standard location, the path to the :term:`PyTrilinos`
    site-packages directory should be added to the :envvar:`PYTHONPATH`
    environment variable; this should be of the form
    :file:`${{INSTALL_DIR}}/lib/${{PYTHON_VERSION}}/site-packages/`. Also,
    the path to the :term:`Trilinos` ``lib`` directory should be added to
    the :envvar:`LD_LIBRARY_PATH` (on Linux) or :envvar:`DYLD_LIBRARY_PATH`
    (on Mac OS X) environment variable; this should be of the form
    :file:`${{INSTALL_DIR}}/lib``.

.. _MPI4PY:

mpi4py
======

https://mpi4py.readthedocs.io/

For :ref:`PARALLEL`, :term:`FiPy` requires ``mpi4py``, in addition to
:term:`PETSc` or :term:`Trilinos`.
