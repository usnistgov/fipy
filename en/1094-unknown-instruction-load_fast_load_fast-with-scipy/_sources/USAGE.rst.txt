.. _USAGE:

==========
Using FiPy
==========

This document explains how to use :term:`FiPy` in a practical sense.
To see the problems that :term:`FiPy` is capable of solving, you can
run any of the scripts in the :ref:`examples <part:examples>`.

.. note::

   We strongly recommend you proceed through the :ref:`examples
   <part:examples>`, but at the very least work through
   :mod:`examples.diffusion.mesh1D` to understand the notation and
   basic concepts of :term:`FiPy`.

We exclusively use either the UNIX command line or :term:`IPython` to
interact with :term:`FiPy`. The commands in the :ref:`examples
<part:examples>` are written with the assumption that they will be
executed from the command line. For instance, from within the main
:term:`FiPy` directory, you can type::

    $ python examples/diffusion/mesh1D.py

A viewer should appear and you should be prompted through a series of
examples.

.. note::

   From within :term:`IPython`, you would type::

       >>> run examples/diffusion/mesh1D.py

In order to customize the examples, or to develop your own scripts, some
knowledge of Python syntax is required.  We recommend you familiarize
yourself with the excellent `Python tutorial`_ :cite:`PythonTutorial`
or with `Dive Into Python`_ :cite:`DiveIntoPython`. Deeper insight into 
Python can be obtained from the :cite:`PythonReference`.

.. _Python tutorial: http://docs.python.org/tut/tut.html
.. _Dive Into Python: http://diveintopython.org

As you gain experience, you may want to browse through the
:ref:`FlagsAndEnvironmentVariables` that affect :term:`FiPy`.

.. _LOGGING:

-------
Logging
-------

Diagnostic information about a :term:`FiPy` run can be obtained using the
:mod:`logging` module.  For example, at the beginning of your script, you
can add::

    >>> import logging
    >>> log = logging.getLogger("fipy")
    >>> console = logging.StreamHandler()
    >>> console.setLevel(logging.INFO)
    >>> log.addHandler(console)

in order to see informational messages in the terminal.  To have more
verbose debugging information save to a file::

    >>> logfile = logging.FileHandler(filename="fipy.log")
    >>> logfile.setLevel(logging.DEBUG)
    >>> log.addHandler(logfile)

    >>> log.setLevel(logging.DEBUG)

To restrict logging to, e.g., information about the :term:`PETSc` solvers::

    >>> petsc = logging.Filter('fipy.solvers.petsc')
    >>> logfile.addFilter(petsc)

More complex configurations can be specified by setting the
:envvar:`FIPY_LOG_CONFIG` environment variable.  In this case, it is not
necessary to add any logging instructions to your own script.  Example
configuration files can be found in
:file:`{FiPySource}/fipy/tools/logging/`.

If `Solving in Parallel`_, the
`mpilogging <https://github.com/usnistgov/mpilogging>`_
package enables reporting which MPI rank each log entry comes from.  For
example::

    >>> from mpilogging import MPIScatteredFileHandler
    >>> mpilog = MPIScatteredFileHandler(filepattern="fipy.%(mpirank)d_of_%(mpisize)d.log"
    >>> mpilog.setLevel(logging.DEBUG)
    >>> log.addHandler(mpilog)

will generate a unique log file for each MPI rank.

------------
Testing FiPy
------------

For a general installation, :term:`FiPy` can be tested by running::

    $ python -c "import fipy; fipy.test()"

This command runs all the test cases in :ref:`FiPy's modules
<part:modules>`, but doesn't include any of the tests in :ref:`FiPy's
examples <part:examples>`. To run the test cases in both :ref:`modules
<part:modules>` and :ref:`examples <part:examples>`, use::

    $ python setup.py test

.. note::

   You may need to first run::

        $ python setup.py egg_info

   for this to work properly.

in an unpacked :term:`FiPy` archive. The test suite can be run with a
number of different configurations depending on which solver suite is
available and other factors. See :ref:`FlagsAndEnvironmentVariables`
for more details.

:term:`FiPy` will skip tests that depend on :ref:`OPTIONALPACKAGES` that
have not been installed. For example, if :term:`Mayavi` and :term:`Gmsh`
are not installed, :term:`FiPy` will warn something like::

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Skipped 131 doctest examples because `gmsh` cannot be found on the $PATH
    Skipped 42 doctest examples because the `tvtk` package cannot be imported
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Although the test suite may show warnings, there should be no other errors.
Any errors should be investigated or reported on the `issue tracker`_.
Users can see if there are any known problems for the latest :term:`FiPy`
distribution by checking :term:`FiPy`'s :ref:`CONTINUOUSINTEGRATION` dashboard.

.. _issue tracker: https://github.com/usnistgov/fipy/issues/new

Below are a number of common `Command-line Flags`_ for testing various
:term:`FiPy` configurations.

Parallel Tests
==============

If :term:`FiPy` is configured for :ref:`PARALLEL`, you can run the tests
on multiple processor cores with::

    $ mpirun -np {# of processors} python setup.py test --trilinos

or::

    $ mpirun -np {# of processors} python -c "import fipy; fipy.test('--trilinos')"

.. _FlagsAndEnvironmentVariables:

--------------------------------------------
Command-line Flags and Environment Variables
--------------------------------------------

:term:`FiPy` chooses a default run time configuration based on the
available packages on the system. The `Command-line Flags`_ and
`Environment Variables`_ sections below describe how to override
:term:`FiPy`'s default behavior.

Command-line Flags
==================

You can add any of the following case-insensitive flags after the name of a
script you call from the command line, e.g.::

    $ python myFiPyScript --someflag

.. cmdoption:: --inline

   Causes many mathematical operations to be performed in C, rather than
   Python, for improved performance. Requires the :mod:`weave`
   package.

.. cmdoption:: --cache

   Causes lazily evaluated :term:`FiPy`
   :class:`~fipy.variables.variable.Variable` objects to retain their
   value.

.. cmdoption:: --no-cache

   Causes lazily evaluated :term:`FiPy`
   :class:`~fipy.variables.variable.Variable` objects to always recalculate
   their value.

The following flags take precedence over the :envvar:`FIPY_SOLVERS`
environment variable:

.. cmdoption:: --pysparse

   Forces the use of the :ref:`PYSPARSE` solvers.

.. cmdoption:: --trilinos

   Forces the use of the :ref:`TRILINOS` solvers, but uses
   :ref:`PYSPARSE` to construct the matrices.

.. cmdoption:: --no-pysparse

   Forces the use of the :ref:`TRILINOS` solvers without any use of
   :ref:`PYSPARSE`.

.. cmdoption:: --scipy

   Forces the use of the :ref:`SCIPY` solvers.

.. cmdoption:: --pyamg

   Forces the use of the :ref:`PYAMG` preconditioners in conjunction
   with the :ref:`SCIPY` solvers.

.. cmdoption:: --pyamgx

   Forces the use of the :ref:`PYAMGX` solvers.

.. cmdoption:: --lsmlib

   Forces the use of the :ref:`LSMLIBDOC` level set solver.

.. cmdoption:: --skfmm

   Forces the use of the :ref:`SCIKITFMM` level set solver.


Environment Variables
=====================

You can set any of the following environment variables in the manner
appropriate for your shell. If you are not running in a shell (*e.g.*,
you are invoking :term:`FiPy` scripts from within :term:`IPython` or IDLE),
you can set these variables via the :data:`os.environ` dictionary,
but you must do so before importing anything from the :mod:`fipy`
package.

.. envvar:: FIPY_DISPLAY_MATRIX

   .. currentmodule:: fipy.terms.term

   If present, causes the graphical display of the solution matrix of each
   equation at each call of :meth:`~Term.solve` or :meth:`~Term.sweep`.
   Setting the value to "``terms``" causes the display of the matrix for each
   :class:`Term` that composes the equation. Requires the :term:`Matplotlib`
   package. Setting the value to "``print``" causes the matrix to be
   printed to the console.

.. envvar:: FIPY_INLINE

   If present, causes many mathematical operations to be performed in C,
   rather than Python. Requires the :mod:`weave` package.

.. envvar:: FIPY_INLINE_COMMENT

   If present, causes the addition of a comment showing the Python context
   that produced a particular piece of :mod:`weave` C code. Useful
   for debugging.

.. envvar:: FIPY_LOG_CONFIG

   Specifies a :term:`JSON`-formatted logging configuration file, suitable
   for passing to :func:`logging.config.dictConfig`.  Example configuration
   files can be found in :file:`{FiPySource}/fipy/tools/logging/`.

.. envvar:: FIPY_SOLVERS

   Forces the use of the specified suite of linear solvers.  Valid
   (case-insensitive) choices are "``petsc``", "``scipy``", "``pysparse``",
   "``trilinos``", "``no-pysparse``", and "``pyamg``".

.. envvar:: FIPY_VERBOSE_SOLVER

   If present, causes the
   :class:`~fipy.solvers.pyAMG.linearGeneralSolver.LinearGeneralSolver` to
   print a variety of diagnostic information.  All other solvers should use
   `Logging`_ and :envvar:`FIPY_LOG_CONFIG`.

.. envvar:: FIPY_VIEWER

   Forces the use of the specified viewer. Valid values are any
   :samp:`{<viewer>}` from the
   :samp:`fipy.viewers.{<viewer>}Viewer`
   modules. The special value of ``dummy`` will allow the script
   to run without displaying anything.

.. envvar:: FIPY_INCLUDE_NUMERIX_ALL

   If present, causes the inclusion of all functions and variables of the
   :mod:`~fipy.tools.numerix` module in the :mod:`fipy` namespace.

.. envvar:: FIPY_CACHE

   If present, causes lazily evaluated :term:`FiPy` 
   :class:`~fipy.variables.variable.Variable` objects to
   retain their value.

.. envvar:: PETSC_OPTIONS

   `PETSc configuration options`_.  Set to "`-help`" and run a script with
   :ref:`PETSC` solvers in order to see what options are possible.  Ignored
   if solver is not :ref:`PETSC`.

.. _PETSc configuration options: https://docs.petsc.org/en/latest/manual/other/#sec-options

.. _PARALLEL:

-------------------
Solving in Parallel
-------------------

:term:`FiPy` can use :ref:`PETSC` or :ref:`TRILINOS` to solve equations in
parallel. Most mesh classes in :mod:`fipy.meshes` can solve in
parallel. This includes all "``...Grid...``" and "``...Gmsh...``"
class meshes. Currently, the only remaining serial-only meshes are
:class:`~fipy.meshes.tri2D.Tri2D` and
:class:`~fipy.meshes.skewedGrid2D.SkewedGrid2D`.

.. attention::

   :term:`FiPy` requires :term:`mpi4py` to work in parallel.

.. tip::

   You are strongly advised to force the use of only one :term:`OpenMP`
   thread with :ref:`Trilinos`::

       $ export OMP_NUM_THREADS=1

   See :ref:`THREADS_VS_RANKS` for more information.

.. note::

   `Trilinos 12.12 has support for Python 3`_, but
   `PyTrilinos on conda-forge`_ presently only provides 12.10, which is
   limited to Python 2.x. :ref:`PETSC` is available for both :term:`Python
   3` and :term:`Python` 2.7.

It should not generally be necessary to change anything in your script.
Simply invoke::

    $ mpirun -np {# of processors} python myScript.py --petsc

or::

    $ mpirun -np {# of processors} python myScript.py --trilinos

instead of::

    $ python myScript.py

The following plot shows the scaling behavior for multiple
processors.  We compare solution time vs number of Slurm_ tasks (available
cores) for a `Method of Manufactured Solutions Allen-Cahn problem`_.

.. plot:: pyplots/scaling.py

   Scaling behavior of different solver packages

"Speedup" relative to :ref:`PySparse` (bigger numbers are better) versus
number of tasks (processes) on a log-log plot.  The number of threads per
:term:`MPI` rank is indicated by the line style (see legend).
:term:`OpenMP` threads :math:`\times` :term:`MPI` ranks = Slurm_ tasks.

A few things can be observed in this plot:

- Both :ref:`PETSc` and :ref:`Trilinos` exhibit power law scaling, but the
  power is only about 0.7.  At least one source of this poor scaling is
  that our "``...Grid...``" meshes parallelize by dividing the mesh into
  slabs, which leads to more communication overhead than more compact
  partitions.  The "``...Gmsh...``" meshes partition more efficiently, but
  carry more overhead in other ways.  We'll be making efforts to improve
  the partitioning of the "``...Grid...``" meshes in a future release.

- :ref:`PETSc` and :ref:`Trilinos` have fairly comparable performance, but
  lag :ref:`PySparse` by a considerable margin.  The :ref:`SciPy` solvers
  are even worse.  Some of this discrepancy may be because the different
  packages are not all doing the same thing.  Different solver packages
  have different default solvers and preconditioners.  Moreover, the
  meaning of the solution tolerance depends on the normalization the solver
  uses and it is not always obvious which of several possibilities a
  particular package employs.  We will endeavor to normalize the
  normalizations in a future release.

- :ref:`PETSc` with one thread is faster than with two threads until the
  number of tasks reaches about 10 and is faster than with four threads
  until the number of tasks reaches more than 20.  :ref:`Trilinos` with one
  thread is faster than with two threads until the number of tasks is more
  than 30.  We don't fully understand the reasons for this, but there may be
  a *modest* benefit, *when using a large number of cpus*, to allow two to
  four :term:`OpenMP` threads per :term:`MPI` rank.  See
  :ref:`THREADS_VS_RANKS` for caveats and more information.

These results are likely both problem and architecture dependent.  You
should develop an understanding of the scaling behavior of your own codes
before doing "production" runs.

The easiest way to confirm that :term:`FiPy` is properly configured to
solve in parallel is to run one of the examples, e.g.,::

    $ mpirun -np 2 examples/diffusion/mesh1D.py

You should see two viewers open with half the simulation running in one of
them and half in the other. If this does not look right (e.g., you get two
viewers, both showing the entire simulation), or if you just want to be
sure, you can run a diagnostic script::

    $ mpirun -np 3 python examples/parallel.py

which should print out::

         mpi4py            PyTrilinos           petsc4py                   FiPy            
    processor 0 of 3 :: processor 0 of 3 :: processor 0 of 3 :: 5 cells on processor 0 of 3
    processor 1 of 3 :: processor 1 of 3 :: processor 1 of 3 :: 7 cells on processor 1 of 3
    processor 2 of 3 :: processor 2 of 3 :: processor 2 of 3 :: 6 cells on processor 2 of 3

If there is a problem with your parallel environment, it should be clear
that there is either a problem importing one of the required packages or
that there is some problem with the :term:`MPI` environment. For example::

         mpi4py            PyTrilinos           petsc4py                   FiPy            
    processor 0 of 3 :: processor 0 of 1 :: processor 0 of 3 :: 10 cells on processor 0 of 1
    [my.machine.com:69815] WARNING: There were 4 Windows created but not freed.
    processor 1 of 3 :: processor 0 of 1 :: processor 1 of 3 :: 10 cells on processor 0 of 1
    [my.machine.com:69814] WARNING: There were 4 Windows created but not freed.
    processor 2 of 3 :: processor 0 of 1 :: processor 2 of 3 :: 10 cells on processor 0 of 1
    [my.machine.com:69813] WARNING: There were 4 Windows created but not freed.

indicates :term:`mpi4py` is properly communicating with :term:`MPI` and is running
in parallel, but that :ref:`TRILINOS` is not, and is running three separate
serial environments. As a result, :term:`FiPy` is limited to three separate
serial operations, too. In this instance, the problem is that although
:ref:`TRILINOS` was compiled with :term:`MPI` enabled, it was compiled against a
different :term:`MPI` library than is currently available (and which :term:`mpi4py`
was compiled against).  The solution, in this instance, is to solve with
:ref:`PETSC` or to rebuild :ref:`TRILINOS` against the active :term:`MPI`
libraries.

When solving in parallel, :term:`FiPy` essentially breaks the problem
up into separate sub-domains and solves them (somewhat) independently.
:term:`FiPy` generally "does the right thing", but if you find that
you need to do something with the entire solution, you can use
``var.``:attr:`~fipy.variables.cellVariable.CellVariable.globalValue`.

.. note:: One option for debugging in parallel is::

   $ mpirun -np {# of processors} xterm -hold -e "python -m ipdb myScript.py"

.. _Trilinos 12.12 has support for Python 3 : https://github.com/trilinos/Trilinos/issues/3203
.. _PyTrilinos on conda-forge: https://anaconda.org/conda-forge/pytrilinos
.. _Slurm: https://slurm.schedmd.com

.. _THREADS_VS_RANKS:

OpenMP Threads vs. MPI Ranks
============================

By default, :ref:`PETSc` and :ref:`Trilinos` spawn as many :term:`OpenMP` threads
as there are cores
available.  This may very well be an intentional optimization, where
they are designed to have one :term:`MPI` rank per node of a cluster, so each
of the child threads would help with computation but would not compete for
I/O resources during ghost cell exchanges and file I/O. However, Python's
`Global Interpreter Lock`_ (GIL) binds all of the child threads to the same
core as their parent!  So instead of improving performance, each core
suffers a heavy overhead from managing those idling threads.

The solution to this is to force these solvers to use only one :term:`OpenMP` thread::

   $ export OMP_NUM_THREADS=1

Because this environment variable affects all processes launched in the
current session, you may prefer to restrict its use to :term:`FiPy` runs::

   $ OMP_NUM_THREADS=1 mpirun -np {# of processors} python myScript.py --trilinos

The difference can be extreme.  We have observed the :term:`FiPy` test
suite to run in `just over two minutes`_ when ``OMP_NUM_THREADS=1``,
compared to `over an hour and 23 minutes`_ when :term:`OpenMP` threads are
unrestricted. We don't know why, but `other platforms`_ do not suffer the
same degree of degradation.

Conceivably, allowing these parallel solvers unfettered access to
:term:`OpenMP` threads with no :term:`MPI` communication at all could perform
as well or better than purely :term:`MPI`
parallelization.  The plot below demonstrates this is not the case.  We
compare solution time vs number of :term:`OpenMP` threads for fixed number of
slots for a `Method of Manufactured Solutions Allen-Cahn problem`_.
:term:`OpenMP` threading always slows down FiPy performance.

.. plot:: pyplots/cpus_vs_threads.py

   Effect of having more :term:`OpenMP` threads for each :term:`MPI` rank

"Speedup" relative to one thread (bigger numbers are better) versus number
of threads for 32 Slurm_ tasks on a log-log plot.  :term:`OpenMP` threads
:math:`\times` :term:`MPI` ranks = Slurm_ tasks.

See https://www.mail-archive.com/fipy@nist.gov/msg03393.html for further
analysis.

It may be possible to configure these packages to use only one :term:`OpenMP` thread,
but this is not the configuration of the version available from conda-forge_
and building Trilinos, at least, is |NotFun (TM)|_.

.. _Global Interpreter Lock:     https://docs.python.org/2.7/c-api/init.html#thread-state-and-the-global-interpreter-lock
.. _just over two minutes:       https://circleci.com/gh/guyer/fipy/461
.. _over an hour and 23 minutes: https://circleci.com/gh/guyer/fipy/423
.. _other platforms:             https://travis-ci.org/usnistgov/fipy/builds/509556033
.. _Method of Manufactured Solutions Allen-Cahn problem:  https://pages.nist.gov/pfhub/benchmarks/benchmark7.ipynb
.. _conda-forge:                 https://conda-forge.github.io/
.. |NotFun (TM)|                 unicode:: NotFun U+2122
.. _NotFun (TM):                 https://commons.wikimedia.org/wiki/File:Hieronymus_Bosch_-_Triptych_of_Garden_of_Earthly_Delights_(detail)_-_WGA2526.jpg#/media/File:Hieronymus_Bosch_-_Triptych_of_Garden_of_Earthly_Delights_(detail)_-_WGA2526.jpg

.. _MeshingWithGmsh:

-----------------
Meshing with Gmsh
-----------------

:term:`FiPy` works with arbitrary polygonal meshes generated by
:term:`Gmsh`.  :term:`FiPy` provides two wrappers classes
(:class:`~fipy.meshes.gmshImport.Gmsh2D` and
:class:`~fipy.meshes.gmshImport.Gmsh3D`) enabling :term:`Gmsh` to be
used directly from python. The classes can be instantiated with a set
of :term:`Gmsh` style commands (see
:mod:`examples.diffusion.circle`). The classes can also be
instantiated with the path to either a :term:`Gmsh` geometry file
(``.geo``) or a :term:`Gmsh` mesh file (``.msh``) (see
:mod:`examples.diffusion.anisotropy`).

As well as meshing arbitrary geometries, :term:`Gmsh` partitions
meshes for parallel simulations. Mesh partitioning automatically
occurs whenever a parallel communicator is passed to the mesh on
instantiation. This is the default setting for all meshes that work in
parallel including :class:`~fipy.meshes.gmshImport.Gmsh2D` and
:class:`~fipy.meshes.gmshImport.Gmsh3D`.

.. note::

    :term:`FiPy` solution accuracy can be compromised with highly
    non-orthogonal or non-conjunctional meshes.

.. _CoupledEquations:

----------------------------
Coupled and Vector Equations
----------------------------

Equations can now be coupled together so that the contributions from
all the equations appear in a single system matrix. This results in
tighter coupling for equations with spatial and temporal derivatives
in more than one variable. In :term:`FiPy` equations are coupled
together using the ``&`` operator::

   >>> eqn0 = ...
   >>> eqn1 = ...
   >>> coupledEqn = eqn0 & eqn1

The ``coupledEqn`` will use a combined system matrix that includes
four quadrants for each of the different variable and equation
combinations. In previous versions of :term:`FiPy` there has been no
need to specify which variable a given term acts on when generating
equations. The variable is simply specified when calling ``solve`` or
``sweep`` and this functionality has been maintained in the case of
single equations. However, for coupled equations the variable that a
given term operates on now needs to be specified when the equation is
generated. The syntax for generating coupled equations has the form::

   >>> eqn0 = Term00(coeff=..., var=var0) + Term01(coeff..., var=var1) == source0
   >>> eqn1 = Term10(coeff=..., var=var0) + Term11(coeff..., var=var1) == source1
   >>> coupledEqn = eqn0 & eqn1

and there is now no need to pass any variables when solving::

   >>> coupledEqn.solve(dt=..., solver=...)

In this case the matrix system will have the form

.. math::

   \left(
   \begin{array}{c|c}
   \text{\ttfamily Term00} & \text{\ttfamily Term01} \\ \hline
   \text{\ttfamily Term10} & \text{\ttfamily Term11}
   \end{array} \right)
   \left(
   \begin{array}{c}
   \text{\ttfamily var0}  \\ \hline
   \text{\ttfamily var1}
   \end{array} \right)
   =
   \left(
   \begin{array}{c}
   \text{\ttfamily source0}  \\ \hline
   \text{\ttfamily source1}
   \end{array} \right)

:term:`FiPy` tries to make sensible decisions regarding each term's
location in the matrix and the ordering of the variable column
array. For example, if ``Term01`` is a transient term then ``Term01``
would appear in the upper left diagonal and the ordering of the
variable column array would be reversed.

The use of coupled equations is described in detail in
:mod:`examples.diffusion.coupled`. Other examples that demonstrate the
use of coupled equations are :mod:`examples.phase.binaryCoupled`,
:mod:`examples.phase.polyxtalCoupled` and
:mod:`examples.cahnHilliard.mesh2DCoupled`. As well as coupling
equations, true vector equations can now be written in :term:`FiPy`.

.. attention::

    Coupled equations are not compatible with
    :ref:`discret-higherOrderDiffusion` terms.  This is not a practical
    limitation, as any higher order terms can be decomposed into multiple
    2nd-order equations.  For example, the pair of `coupled Cahn-Hilliard &
    Allen-Cahn`_ 4th- and 2nd-order equations

    .. math::

        \frac{\partial C}{\partial t}
        &= \nabla\cdot\left[
            M\nabla\left(
                \frac{\partial f(c, \phi)}{\partial C}
                - \kappa_C\nabla^2 C
            \right)
        \right]
        \\
        \frac{\partial \phi}{\partial t}
        &= -L\left(
            \frac{\partial f(c, \phi)}{\partial \phi}
            - \kappa_\phi\nabla^2 \phi
        \right)

    can be decomposed to three 2nd-order equations

    .. math::

        \frac{\partial C}{\partial t}
        &= \nabla\cdot\left(
            M\nabla\mu
        \right)
        \\
        \mu
        &= \frac{\partial f(c, \phi)}{\partial C}
           - \kappa_C\nabla^2 C
        \\
        \frac{\partial \phi}{\partial t}
        &= -L\left(
            \frac{\partial f(c, \phi)}{\partial \phi}
            - \kappa_\phi\nabla^2 \phi
        \right)

    .. _coupled Cahn-Hilliard & Allen-Cahn: https://pages.nist.gov/pfhub/benchmarks/benchmark2.ipynb

.. _BoundaryConditions:

-------------------
Boundary Conditions
-------------------

.. currentmodule:: fipy.variables.cellVariable

Default boundary conditions
===========================

If no constraints are applied, solutions are conservative, i.e., all
boundaries are zero flux.  For the equation

.. math::

    \frac{\partial\phi}{\partial t}
    &= \nabla\cdot\left(\vec{a}\phi\right) + \nabla\cdot\left(b\nabla\phi\right)

the condition on the boundary :math:`S` is

.. math::

   \hat{n}\cdot\left(\vec{a}\phi + b\nabla\phi\right) = 0\qquad\text{on $S$.}

Applying fixed value (Dirichlet) boundary conditions
====================================================

To apply a fixed value boundary condition use the
:meth:`~CellVariable.constrain` method. For example, to fix `var` to
have a value of `2` along the upper surface of a domain, use

>>> var.constrain(2., where=mesh.facesTop)

.. note::

   The old equivalent
   :class:`~fipy.boundaryConditions.fixedValue.FixedValue` boundary
   condition is now deprecated.

Applying fixed gradient boundary conditions (Neumann)
=====================================================

To apply a fixed Gradient boundary condition use the
:attr:`~.CellVariable.faceGrad`.\
:meth:`~fipy.variables.variable.Variable.constrain` method. For
example, to fix `var` to have a gradient of `(0,2)` along the upper
surface of a 2D domain, use

>>> var.faceGrad.constrain(((0,),(2,)), where=mesh.facesTop)

If the gradient normal to the boundary (*e.g.*,
:math:`\hat{n}\cdot\nabla\phi`) is to be set to a scalar value of `2`, use

>>> var.faceGrad.constrain(2 * mesh.faceNormals, where=mesh.exteriorFaces)

Applying fixed flux boundary conditions
=======================================

Generally these can be implemented with a judicious use of
:attr:`~.CellVariable.faceGrad`.\
:meth:`~fipy.variables.variable.Variable.constrain`.  Failing that, an
exterior flux term can be added to the equation. Firstly, set the
terms' coefficients to be zero on the exterior faces,

>>> diffCoeff.constrain(0., mesh.exteriorFaces)
>>> convCoeff.constrain(0., mesh.exteriorFaces)

then create an equation with an extra term to account for the exterior flux,

>>> eqn = (TransientTerm() + ConvectionTerm(convCoeff)
...        == DiffusionCoeff(diffCoeff)
...        + (mesh.exteriorFaces * exteriorFlux).divergence)

where `exteriorFlux` is a rank 1
:class:`~fipy.variables.faceVariable.FaceVariable`.

.. note::

   The old equivalent :class:`~fipy.boundaryConditions.fixedFlux.FixedFlux`
   boundary condition is now deprecated.

Applying outlet or inlet boundary conditions
============================================

Convection terms default to a no flux boundary condition unless the
exterior faces are associated with a constraint, in which case either
an inlet or an outlet boundary condition is applied depending on the
flow direction.

Applying spatially varying boundary conditions
==============================================

The use of spatial varying boundary conditions is best demonstrated with an
example. Given a 2D equation in the domain :math:`0 < x < 1` and :math:`0 < y < 1` with
boundary conditions,

.. math::

  \phi = \left\{
            \begin{aligned}
                xy &\quad \text{on $x>1/2$ and $y>1/2$} \\
                \vec{n} \cdot \vec{F} = 0 &\quad \text{elsewhere}
            \end{aligned}
        \right.

where :math:`\vec{F}` represents the flux. The boundary conditions in :term:`FiPy` can
be written with the following code,

>>> X, Y = mesh.faceCenters
>>> mask =  ((X < 0.5) | (Y < 0.5))
>>> var.faceGrad.constrain(0, where=mesh.exteriorFaces & mask)
>>> var.constrain(X * Y, where=mesh.exteriorFaces & ~mask)

then

>>> eqn.solve(...)

Further demonstrations of spatially varying boundary condition can be found
in :mod:`examples.diffusion.mesh20x20`
and :mod:`examples.diffusion.circle`

Applying Robin boundary conditions
==================================

The Robin condition applied on the portion of the boundary :math:`S_R`

.. math::

   \hat{n}\cdot\left(\vec{a}\phi + b\nabla\phi\right) = g\qquad\text{on $S_R$}

can often be substituted for the flux in an equation

.. math::

    \frac{\partial\phi}{\partial t}
    &= \nabla\cdot\left(\vec{a}\phi\right) + \nabla\cdot\left(b\nabla\phi\right)
    \\
    \int_V\frac{\partial\phi}{\partial t}\,dV
    &= \int_S \hat{n} \cdot \left(\vec{a}\phi + b\nabla\phi\right) \, dS
    \\
    \int_V\frac{\partial\phi}{\partial t}\,dV
    &= \int_{S \notin S_R} \hat{n} \cdot \left(\vec{a}\phi + b\nabla\phi\right) \, dS
    + \int_{S \in S_R} g \, dS

At faces identified by ``mask``,

>>> a = FaceVariable(mesh=mesh, value=..., rank=1)
>>> a.setValue(0., where=mask)
>>> b = FaceVariable(mesh=mesh, value=..., rank=0)
>>> b.setValue(0., where=mask)
>>> g = FaceVariable(mesh=mesh, value=..., rank=0)
>>> eqn = (TransientTerm() == PowerLawConvectionTerm(coeff=a)
...        + DiffusionTerm(coeff=b)
...        + (g * mask * mesh.faceNormals).divergence)

When the Robin condition does not exactly map onto the boundary flux, we
can attempt to apply it term by term.  The Robin condition relates the
gradient at a boundary face to the value on that face, however
:term:`FiPy` naturally calculates variable values at cell centers
and gradients at intervening faces. Using a first order upwind
approximation, the boundary value of the variable at face :math:`f` can be written in terms of
the value at the neighboring cell :math:`P` and the normal gradient at the boundary:

.. math::
   :label: upwind1

   \phi_f &\approx \phi_P + \left(\vec{d}_{Pf}\cdot\nabla\phi\right)_f
   \\
   &\approx \phi_P + \left(\hat{n}\cdot\nabla\phi\right)_f\left(\vec{d}_{Pf}\cdot\hat{n}\right)_f

where :math:`\vec{d}_{Pf}` is the distance vector to the center of the face
:math:`f` from the center of the adjoining cell :math:`P`.  The
approximation :math:`\left(\vec{d}_{Pf}\cdot\nabla\phi\right)_f \approx
\left(\hat{n}\cdot\nabla\phi\right)_f\left(\vec{d}_{Pf}\cdot\hat{n}\right)_f`
is most valid when the mesh is orthogonal.

Substituting this expression into the Robin condition:

.. math::
   :label: Robin_facegrad

   \hat{n}\cdot\left(\vec{a} \phi + b \nabla\phi\right)_f &= g \\
   \hat{n}\cdot\left[\vec{a} \phi_P
   + \vec{a} \left(\hat{n}\cdot\nabla\phi\right)_f\left(\vec{d}_{Pf}\cdot\hat{n}\right)_f
   + b \nabla\phi\right]_f &\approx g \\
   \left(\hat{n}\cdot\nabla\phi\right)_f
   &\approx \frac{g_f - \left(\hat{n}\cdot\vec{a}\right)_f \phi_P}
                 {\left(\vec{d}_{Pf}\cdot\vec{a}\right)_f + b_f}

we obtain an expression for the gradient at the boundary face in terms of
its neighboring cell.  We can, in turn, substitute this back into
:eq:`upwind1`

.. math::
   :label: upwind2

   \phi_f &\approx \phi_P
   + \frac{g_f - \left(\hat{n}\cdot\vec{a}\right)_f \phi_P}
          {\left(\vec{d}_{Pf}\cdot\vec{a}\right)_f + b_f}
   \left(\vec{d}_{Pf}\cdot\hat{n}\right)_f \\
   &\approx \frac{g_f \left(\hat{n}\cdot\vec{d}_{Pf}\right)_f + b_f\phi_P}
                 {\left(\vec{d}_{Pf}\cdot\vec{a}\right)_f + b_f}

to obtain the value on the boundary face in terms of the neighboring cell.

Substituting :eq:`Robin_facegrad` into the discretization of the
:class:`~fipy.terms.diffusionTerm.DiffusionTerm`:

.. math::

   \int_V \nabla\cdot\left(\Gamma\nabla\phi\right) dV
   &= \int_S \Gamma \hat{n}\cdot\nabla\phi\, S \\
   &\approx \sum_f \Gamma_f \left(\hat{n}\cdot\nabla\phi\right)_f A_f \\
   &= \sum_{f \notin S_R} \Gamma_f \left(\hat{n}\cdot\nabla\phi\right)_f A_f
   + \sum_{f \in S_R} \Gamma_f \left(\hat{n}\cdot\nabla\phi\right)_f A_f \\
   &\approx \sum_{f \notin S_R} \Gamma_f \left(\hat{n}\cdot\nabla\phi\right)_f A_f
   + \sum_{f \in S_R} \Gamma_f \frac{g_f - \left(\hat{n}\cdot\vec{a}\right)_f \phi_P}
                       {\left(\vec{d}_{Pf}\cdot\vec{a}\right)_f + b_f} A_f

An equation of the form

>>> eqn = TransientTerm() == DiffusionTerm(coeff=Gamma0)

can be constrained to have a Robin condition at faces identified by
``mask`` by making the following modifications

>>> Gamma = FaceVariable(mesh=mesh, value=Gamma0)
>>> Gamma.setValue(0., where=mask)
>>> dPf = FaceVariable(mesh=mesh, 
...                    value=mesh._faceToCellDistanceRatio * mesh.cellDistanceVectors)
>>> n = mesh.faceNormals
>>> a = FaceVariable(mesh=mesh, value=..., rank=1)
>>> b = FaceVariable(mesh=mesh, value=..., rank=0)
>>> g = FaceVariable(mesh=mesh, value=..., rank=0)
>>> RobinCoeff = mask * Gamma0 * n / (dPf.dot(a) + b)
>>> eqn = (TransientTerm() == DiffusionTerm(coeff=Gamma) + (RobinCoeff * g).divergence
...        - ImplicitSourceTerm(coeff=(RobinCoeff * n.dot(a)).divergence)

Similarly, for a :class:`~fipy.terms.convectionTerm.ConvectionTerm`, we can
substitute :eq:`upwind2`:

.. math::

   \int_V \nabla\cdot\left(\vec{u}\phi\right) dV
   &= \int_S \hat{n}\cdot\vec{u} \phi\,dS \\
   &\approx \sum_f \left(\hat{n}\cdot\vec{u}\right)_f \phi_f A_f \\
   &= \sum_{f \notin S_R} \left(\hat{n}\cdot\vec{u}\right)_f \phi_f A_f
   + \sum_{f \in S_R} \left(\hat{n}\cdot\vec{u}\right)_f
        \frac{g_f \left(\hat{n}\cdot\vec{d}_{Pf}\right)_f + b_f\phi_P}
             {\left(\vec{d}_{Pf}\cdot\vec{a}\right)_f + b_f} A_f

.. note:: An expression like the heat flux convection boundary condition
   :math:`-k\nabla T\cdot\hat{n} = h(T - T_\infty)` can be put in the form of the
   Robin condition used above by letting :math:`\vec{a} \equiv h \hat{n}`,
   :math:`b \equiv k`, and :math:`g \equiv h T_\infty`.

Applying internal "boundary" conditions
=======================================

Applying internal boundary conditions can be achieved through the use
of implicit and explicit sources. 

Internal fixed value
--------------------

An equation of the form

>>> eqn = TransientTerm() == DiffusionTerm()

can be constrained to have a fixed internal ``value`` at a position
given by ``mask`` with the following alterations

>>> eqn = (TransientTerm() == DiffusionTerm() 
...                           - ImplicitSourceTerm(mask * largeValue) 
...                           + mask * largeValue * value)

The parameter ``largeValue`` must be chosen to be large enough to
completely dominate the matrix diagonal and the RHS vector in cells
that are masked. The ``mask`` variable would typically be a
``CellVariable`` Boolean constructed using the cell center values.

Internal fixed gradient
-----------------------

An equation of the form

>>> eqn = TransientTerm() == DiffusionTerm(coeff=Gamma0)

can be constrained to have a fixed internal ``gradient`` magnitude 
at a position given by ``mask`` with the following alterations

>>> Gamma = FaceVariable(mesh=mesh, value=Gamma0)
>>> Gamma[mask.value] = 0.
>>> eqn = (TransientTerm() == DiffusionTerm(coeff=Gamma) 
...        + DiffusionTerm(coeff=largeValue * mask)
...        - ImplicitSourceTerm(mask * largeValue * gradient 
...                             * mesh.faceNormals).divergence)

The parameter ``largeValue`` must be chosen to be large enough to
completely dominate the matrix diagonal and the RHS vector in cells
that are masked. The ``mask`` variable would typically be a
``FaceVariable`` Boolean constructed using the face center values.

Internal Robin condition
------------------------

Nothing different needs to be done when 
`applying Robin boundary conditions`_ at internal faces.

.. note::

   While we believe the derivations for
   `applying Robin boundary conditions`_ are "correct", they often do not
   seem to produce the intuitive result. At this point, we think this has 
   to do with the pathology of "internal" boundary conditions, but remain 
   open to other explanations. :term:`FiPy` was designed with diffuse 
   interface treatments (phase field and level set) in mind and, as such, 
   internal "boundaries" do not come up in our own work and have not 
   received much attention.

.. warning::

  The constraints mechanism is not designed to constrain internal values
  for variables that are being solved by equations. In particular, one must
  be careful to distinguish between constraining internal cell values
  during the solve step and simply applying arbitrary constraints to a
  ``CellVariable``. Applying a constraint,

  >>> var.constrain(value, where=mask)

  simply fixes the returned value of ``var`` at ``mask`` to be
  ``value``. It does not have any effect on the implicit value of ``var`` at the
  ``mask`` location during the linear solve so it is not a substitute
  for the source term machinations described above. Future releases of
  :term:`FiPy` may implicitly deal with this discrepancy, but the current
  release does not. 

  A simple example can be used to demonstrate this::

  >>> m = Grid1D(nx=2, dx=1.)
  >>> var = CellVariable(mesh=m)

  We wish to solve :math:`\nabla^2 \phi = 0` subject to
  :math:`\phi\rvert_\text{right} = 1` and :math:`\phi\rvert_{x < 1} = 0.25`.
  We apply a constraint to the faces for the right side boundary condition
  (which works).

  >>> var.constrain(1., where=m.facesRight)

  We create the equation with the source term constraint described above

  >>> mask = m.x < 1.
  >>> largeValue = 1e+10
  >>> value = 0.25
  >>> eqn = DiffusionTerm() - ImplicitSourceTerm(largeValue * mask) + largeValue * mask * value

  and the expected value is obtained.

  >>> eqn.solve(var)
  >>> print var
  [ 0.25  0.75]

  However, if a constraint is used without the source term constraint an
  unexpected solution is obtained

  >>> var.constrain(0.25, where=mask)
  >>> eqn = DiffusionTerm()
  >>> eqn.solve(var)
  >>> print var
  [ 0.25  1.  ]

  although the left cell has the expected value as it is constrained.

  :term:`FiPy` has simply solved :math:`\nabla^2 \phi = 0` with
  :math:`\phi\rvert_\text{right} = 1` and (by default)
  :math:`\hat{n}\cdot\nabla\phi\rvert_\text{left} = 0`, giving :math:`\phi
  = 1` everywhere, and then subsequently replaced the cells :math:`x < 1`
  with :math:`\phi = 0.25`.

.. %    http://thread.gmane.org/gmane.comp.python.fipy/726
   %    http://thread.gmane.org/gmane.comp.python.fipy/846

   %    \subsection{Fourth order boundary conditions}

   %    http://thread.gmane.org/gmane.comp.python.fipy/923

   %    \subsection{Periodic boundary conditions}

   %    http://thread.gmane.org/gmane.comp.python.fipy/135

   %    \subsection{Time dependent boundary conditions}

   %    http://thread.gmane.org/gmane.comp.python.fipy/2

   %    \subsection{Internal boundary conditions}

.. _RunningUnderPython2:

----------------------
Running under Python 2
----------------------

Thanks to the future_ package and to the contributions of pya_ and
woodscn_, :term:`FiPy` runs under both :term:`Python 3` and :term:`Python`
2.7, without conversion or modification.

Because :term:`Python` itself will `drop support for Python 2.7 on January
1, 2020`_ and many of the prerequisites for :term:`FiPy` have `pledged to
drop support for Python 2.7 no later than 2020`_, we have prioritized adding
support for better :term:`Python 3` solvers, starting with
:term:`petsc4py`.

Because the faster :term:`PySparse` and :term:`Trilinos` solvers are not
available under :term:`Python 3`, we will maintain :term:`Python` 2.x
support as long as practical. Be aware that the conda-forge_ packages that
:term:`FiPy` depends upon are not well-maintained on :term:`Python` 2.x
and our support for that configuration is rapidly becoming impractical,
despite the present performance benefits. Hopefully, we will learn
how to optimize our use of :ref:`PETSc` and/or :ref:`Trilinos` 12.12 will
become available on conda-forge_.

.. _future: http://python-future.org
.. _pya: https://github.com/pya
.. _woodscn: https://github.com/pya
.. _drop support for Python 2.7 on January 1, 2020: https://www.python.org/dev/peps/pep-0373/#update
.. _pledged to drop support for Python 2.7 no later than 2020: https://python3statement.org

------
Manual
------

You can view the manual online at <http://pages.nist.gov/fipy>.
Alternatively, it may be possible to build a fresh copy by issuing the
following command in the :file:`docs/` directory::

    $ make html

or::

    $ make latexpdf

.. note::

   This mechanism is intended primarily for the developers.  At a minimum,
   you will need `Sphinx <http://www.sphinx-doc.org/>`_, plus all of its
   prerequisites.  We are currently building with Sphinx v7.0.  Python 2.7
   probably won't work.

   We install via conda::

   $ conda install --channel conda-forge sphinx

   Bibliographic citations require the `sphinxcontrib-bibtex` package::

   $ python -m pip install sphinxcontrib-bibtex

   Some documentation uses `numpydoc` styling::

   $ python -m pip install numpydoc

   Some embeded figures require `matplotlib`, `pandas`, and `imagemagick`::

   $ conda install --channel conda-forge matplotlib pandas imagemagick

   The PDF file requires `SIunits.sty`_ available, e.g., from
   `texlive-science`_.

   Spelling is checked automatically in the course of
   :ref:`CONTINUOUSINTEGRATION`. If you wish to check manually, you will
   need `pyspelling`, `hunspell`, and the `libreoffice` dictionaries::

   $ conda install --channel conda-forge hunspell
   $ python -m pip install pyspelling
   $ wget -O en_US.aff  https://cgit.freedesktop.org/libreoffice/dictionaries/plain/en/en_US.aff?id=a4473e06b56bfe35187e302754f6baaa8d75e54f
   $ wget -O en_US.dic https://cgit.freedesktop.org/libreoffice/dictionaries/plain/en/en_US.dic?id=a4473e06b56bfe35187e302754f6baaa8d75e54f

.. _download the latest manual:  http://www.ctcms.nist.gov/fipy/download/
.. _SIunits.sty: https://ctan.org/pkg/siunits
.. _texlive-science: https://packages.debian.org/stretch/texlive-science
