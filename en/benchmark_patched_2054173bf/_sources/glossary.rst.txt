.. _glossary:

Glossary
--------

.. glossary::

   AppVeyor
      A cloud-based :term:`Continuous Integration` tool.
      See https://www.appveyor.com.

   Azure
      A cloud-based :term:`Continuous Integration` tool.
      See https://dev.azure.com.

   Buildbot
      The Buildbot is a system to automate the compile/test cycle
      required by most software projects to validate code changes.
      No longer used for :term:`FiPy`.
      See http://trac.buildbot.net/.

   CircleCI
      A cloud-based :term:`Continuous Integration` tool.
      See https://circleci.com.

   conda
      An open source package management system and environment management
      system that runs on Windows, macOS and Linux.  Conda quickly
      installs, runs and updates packages and their dependencies.  Conda
      easily creates, saves, loads and switches between environments on
      your local computer.  It was created for Python programs, but it can
      package and distribute software for any language.
      See https://conda.io.

   Continuous Integration
      The practice of frequently testing and integrating one's new or
      changed code with the existing code repository.
      See https://en.wikipedia.org/wiki/Continuous_integration.

   FiPy
      The eponymous software package.
      See http://www.ctcms.nist.gov/fipy.

   GitHub Actions
      A cloud-based :term:`Continuous Integration` tool.
      See https://github.com/features/actions.

   Gmsh
      A free and Open Source 3D (and 2D!) finite element grid generator. It
      also has a CAD engine and post-processor that :term:`FiPy` does not
      make use of. See http://www.geuz.org/gmsh.

   IPython
      An improved :term:`Python` shell that integrates nicely with
      :ref:`MATPLOTLIB`. See http://ipython.scipy.org/.

   JSON
      JavaScript Object Notation.  A text format suitable for storing
      structured information such as :class:`dict` or :class:`list`.
      See https://www.json.org/.

   linux
      An operating system.
      See http://www.linux.org.

   macOS
      An operating system.
      See http://www.apple.com/macos.

   Matplotlib
      :mod:`matplotlib` :term:`Python` package displays publication quality
      results. It displays both 1D X-Y type plots and 2D contour plots for
      structured data. It does not display unstructured 2D data or 3D data.
      It works on all common platforms and produces publication quality hard
      copies. See
      http://matplotlib.sourceforge.net
      and :ref:`Matplotlib`.

   Mayavi
      The :mod:`mayavi` Data Visualizer is a free, easy to use scientific data
      visualizer.  It displays 1D, 2D and 3D data. It is the only :term:`FiPy`
      viewer available for 3D data. Other viewers are probably better for 1D
      or 2D viewing. See
      http://code.enthought.com/projects/mayavi
      and :ref:`MAYAVI`.

   MayaVi
      The predecessor to :term:`Mayavi`. Yes, it's confusing.

   MPI
      The Message Passing Interface is a standard that allows the use of
      multiple processors. See
      http://www.mpi-forum.org

   mpi4py
      MPI for Python provides bindings of the Message Passing Interface
      (:term:`MPI`) standard for the Python programming language, allowing
      any Python program to exploit multiple processors.  For
      :ref:`PARALLEL`, :term:`FiPy` requires this package, in addition to
      :term:`PETSc` or :term:`Trilinos`.  See
      https://mpi4py.readthedocs.io.

   numarray
      An archaic predecessor to :term:`NumPy`.

   Numeric
      An archaic predecessor to :term:`NumPy`.

   NumPy
      The :mod:`numpy` :term:`Python` package provides array arithmetic
      facilities. See
      http://www.scipy.org/NumPy.

   OpenMP
      The Open Multi-Processing architecture is a specification for a set
      of compiler directives, library routines, and environment variables
      that can be used to specify high-level parallelism in Fortran and
      C/C++ programs. See
      https://www.openmp.org.

   pandas
      "Python Data Analysis Library" provides high-performance data structures
      for flexible, extensible analysis. See http://pandas.pydata.org.

   PETSc
      The Portable, Extensible Toolkit for Scientific Computation is a
      suite of data structures and routines for the scalable (parallel)
      solution of scientific applications modeled by partial differential
      equations. See https://www.mcs.anl.gov/petsc and :ref:`PETSC`.

   petsc4py
      :term:`Python` wrapper for :term:`PETSc`. See
      https://petsc4py.readthedocs.io/.

   pip
      "pip installs python" is a tool for installing and managing Python
      packages, such as those found in :term:`PyPI`.
      See http://www.pip-installer.org.

   PyAMG
      A suite of python-based preconditioners. See
      http://code.google.com/p/pyamg/
      and :ref:`PYAMG`.

   pyamgx
      a :term:`Python` interface to the NVIDIA 
      `AMGX <https://github.com/NVIDIA/AMGX>`_ library, which can be used
      to construct complex solvers and preconditioners to solve sparse
      sparse linear systems on the GPU. See https://pyamgx.readthedocs.io/
      and :ref:`PYAMGX`.

   PyPI
      The Python Package Index is a repository of software for the
      :term:`Python` programming language.
      See http://pypi.python.org/pypi.

   Pyrex
      A mechanism for mixing C and Python code.
      See http://www.cosc.canterbury.ac.nz/greg.ewing/python/Pyrex/.

   Pysparse
      The :mod:`pysparse` :term:`Python` package provides sparse
      matrix storage, solvers, and linear algebra routines. See
      http://pysparse.sourceforge.net
      and :ref:`PYSPARSE`.

   Python
      The programming language that :term:`FiPy` (and your scripts) are
      written in. See
      http://www.python.org/.

   Python 3
      The (likely) future of the :term:`Python` programming language.
      Third-party packages are slowly being adapted, but many that
      :term:`FiPy` uses are not yet available. See
      http://docs.python.org/py3k/
      and :pep:`3000`.

   PyTrilinos
      :term:`Python` wrapper for :term:`Trilinos`. See
      http://trilinos.sandia.gov/packages/pytrilinos/.

   PyxViewer
      A now defunct python viewer.

   ScientificPython
      A collection of useful utilities for scientists. See
      http://dirac.cnrs-orleans.fr/plone/software/scientificpython.

   SciPy
      The :mod:`scipy` package provides a wide range of scientific and
      mathematical operations. :term:`FiPy` can use
      :term:`Scipy`'s solver suite for linear solutions. See
      http://www.scipy.org/.
      and :ref:`SCIPY`.

   Sphinx
      The tools used to generate the :term:`FiPy` documentation.
      See
      http://sphinx.pocoo.org/.

   TravisCI
      A cloud-based :term:`Continuous Integration` tool.
      See https://travis-ci.org.

   Trilinos
      This package provides sparse matrix storage, solvers, and
      preconditioners, and can be used instead of :term:`Pysparse`.
      :term:`Trilinos` preconditioning allows for iterative solutions
      to some difficult problems that :term:`Pysparse` cannot
      solve. See
      http://trilinos.sandia.gov
      and :ref:`TRILINOS`.

   Weave
      The :mod:`weave` package can enhance performance with C language
      inlining.  See https://github.com/scipy/weave.

   Windows
      An operating system.
      See http://www.microsoft.com/windows.

