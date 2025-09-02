.. _glossary:

Glossary
--------

.. glossary::

   AMG
      Algebraic multigrid method for solving sparse matrices.
      AMG solves linear systems with a general non-symmetric coefficient
      matrix.
      See https://en.wikipedia.org/wiki/Multigrid_method#Algebraic_multigrid_(AMG).

   AppVeyor
      A cloud-based :term:`Continuous Integration` tool.
      See https://www.appveyor.com.

   Azure
      A cloud-based :term:`Continuous Integration` tool.
      See https://dev.azure.com.

   BiCG
      Biconjugate gradient method for solving sparse matrices.
      BiCG solves linear systems with a general non-symmetric coefficient
      matrix.
      See https://en.wikipedia.org/wiki/Biconjugate_gradient_method.

   BiCGSTAB
      Biconjugate gradient (stabilized) method for solving sparse matrices.
      BiCGSTAB solves linear systems with a general non-symmetric
      coefficient matrix.
      See https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method.

   CG
      Conjugate Gradient method for solving sparse matrices.  CG solves
      linear systems with a symmetric positive-definite coefficient matrix.
      See https://en.wikipedia.org/wiki/Conjugate_gradient_method.

   CGS
      Conjugate Gradient Squared method for solving sparse matrices, a
      variant of :term:`BiCG`.  CGS solves linear systems with a general
      non-symmetric coefficient matrix.
      See https://en.wikipedia.org/wiki/Conjugate_gradient_squared_method.

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

   FGMRES
      Flexible Inner-Outer Preconditioned :term:`GMRES` algorithm for
      solving sparse matrices.  FGMRES solves systems with a general
      non-symmetric coefficient matrix.
      See https://doi.org/10.1137/0914028.

   GitHub Actions
      A cloud-based :term:`Continuous Integration` tool.
      See https://github.com/features/actions.

   GMRES
      Generalized Minimal RESidual method for solving sparse matrices.
      GMRES solves systems with a general non-symmetric coefficient matrix.
      See https://en.wikipedia.org/wiki/Generalized_minimal_residual_method.

   Gmsh
      A free and Open Source 3D (and 2D!) finite element grid generator. It
      also has a CAD engine and post-processor that :term:`FiPy` does not
      make use of. See http://www.geuz.org/gmsh.

   IPython
      An improved :term:`Python` shell that integrates nicely with
      :ref:`MATPLOTLIB`. See http://ipython.scipy.org/.

   JOR
      Jacobi over-relaxation method for solving sparse matrices.  JOR
      solves systems with a general non-symmetric coefficient matrix.

   JSON
      JavaScript Object Notation.  A text format suitable for storing
      structured information such as :class:`dict` or :class:`list`.
      See https://www.json.org/.

   linux
      An operating system.
      See http://www.linux.org.

   LU
      Lower-Upper decomposition method for solving sparse matrices.  LU
      solves systems with a general non-symmetric coefficient matrix using
      partial pivoting.
      See https://en.wikipedia.org/wiki/LU_decomposition.

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
      :ref:`PETSC` or :ref:`TRILINOS`.  See
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

   PCG
      Preconditioned conjugate gradient method for solving sparse matrices.
      PCG solves systems with a symmetric positive definite coefficient matrix.
      See https://en.wikipedia.org/wiki/Conjugate_gradient_method.

   petsc4py
      :term:`Python` wrapper for :ref:`PETSC`. See
      https://petsc4py.readthedocs.io/.

   pip
      "pip installs python" is a tool for installing and managing Python
      packages, such as those found in :term:`PyPI`.
      See http://www.pip-installer.org.

   PyPI
      The Python Package Index is a repository of software for the
      :term:`Python` programming language.
      See http://pypi.python.org/pypi.

   Pyrex
      A mechanism for mixing C and Python code.
      See http://www.cosc.canterbury.ac.nz/greg.ewing/python/Pyrex/.

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
      :term:`Python` wrapper for :ref:`TRILINOS`. See
      https://trilinos.github.io/pytrilinos.html.

   PyxViewer
      A now defunct python viewer.

   ScientificPython
      A collection of useful utilities for scientists. See
      http://dirac.cnrs-orleans.fr/plone/software/scientificpython.

   SciPy
      The :mod:`scipy` package provides a wide range of scientific and
      mathematical operations. :term:`FiPy` can use
      :ref:`SCIPY` solver suite for linear solutions. See
      http://www.scipy.org/.

   Sphinx
      The tools used to generate the :term:`FiPy` documentation.
      See
      http://sphinx.pocoo.org/.

   steppyngstounes
      This package provides iterators that simplify both deterministic and
      adaptive stepping in time (or other independent variables).
      See https://pages.nist.gov/steppyngstounes/en/latest.

   TravisCI
      A cloud-based :term:`Continuous Integration` tool.
      See https://travis-ci.org.

   Weave
      The :mod:`weave` package can enhance performance with C language
      inlining.  See https://github.com/scipy/weave.

   Windows
      An operating system.
      See http://www.microsoft.com/windows.

