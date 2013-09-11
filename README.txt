========
Overview
========

.. only:: latex

   :term:`FiPy` is an object oriented, partial differential equation (PDE)
   solver, written in :term:`Python`, based on a standard finite volume
   (FV) approach. The framework has been developed in the Materials Science
   and Engineering Division (MSED_) and Center for Theoretical and
   Computational Materials Science (CTCMS_), in the Material Measurement
   Laboratory (MML_) at the National Institute of Standards and Technology
   (NIST_).

   The solution of coupled sets of PDEs is ubiquitous to the numerical
   simulation of science problems.  Numerous PDE solvers exist, using a
   variety of languages and numerical approaches. Many are proprietary,
   expensive and difficult to customize.  As a result, scientists spend
   considerable resources repeatedly developing limited tools for
   specific problems.  Our approach, combining the FV method and :term:`Python`,
   provides a tool that is extensible, powerful and freely available. A
   significant advantage to :term:`Python` is the existing suite of tools for
   array calculations, sparse matrices and data rendering. 

The :term:`FiPy` framework includes terms for transient diffusion,
convection and standard sources, enabling the solution of arbitrary
combinations of coupled elliptic, hyperbolic and parabolic PDEs. Currently
implemented models include phase field :cite:`BoettingerReview:2002`
:cite:`ChenReview:2002` :cite:`McFaddenReview:2002` treatments of polycrystalline,
dendritic, and electrochemical phase transformations, as well as drug
eluting stents :cite:`Saylor:2011p2794`, reactive wetting :cite:`PhysRevE.82.051601`,
photovoltaics :cite:`Hangarter:2011p2795` and a level set treatment of the
electrodeposition process :cite:`NIST:damascene:2001`.

.. only:: latex
  
   The latest information about :term:`FiPy` can be found at
   http://www.ctcms.nist.gov/fipy/.

---------------------------------
Even if you don't read manuals...
---------------------------------

...please read :ref:`INSTALLATION`, :ref:`USAGE` and :ref:`FAQ`, as well 
as :mod:`examples.diffusion.mesh1D`.

--------------------------------
What's new in version |release|?
--------------------------------

The bump in major version number reflects more on the substantial increase
in capabilities and ease of use than it does on a break in compatibility
with FiPy 2.x. Few, if any, changes to your existing scripts should be
necessary.

The significant changes since version 2.1 are:

- :ref:`CoupledEquations` are now supported.
- A more robust mechanism for specifying :ref:`BoundaryConditions` is now 
  used.
- Most :class:`~fipy.meshes.mesh.Mesh`\es can be partitioned by 
  :ref:`MeshingWithGmsh`.
- :ref:`PYAMG` and :ref:`SCIPY` have been added to the :ref:`SOLVERS`.
- FiPy is capable of :ref:`RunningUnderPython3`.
- "getter" and "setter" methods have been pervasively changed to Python 
  properties.
- The test suite now runs much faster.
- Tests can now be run on a full install using `fipy.test()`.
- The functions of the :mod:`~fipy.tools.numerix` module are no longer 
  included in the :mod:`fipy` namespace. See :mod:`examples.updating.update2_0to3_0` 
  for details.
- Equations containing a :class:`~fipy.terms.transientTerm.TransientTerm`,
  must specify the timestep by passing a ``dt=`` argument when calling
  :meth:`~fipy.terms.term.Term.solve` or :meth:`~fipy.terms.term.Term.sweep`.

Tickets fixed in this release::

    45  Navier Stokes
    85  CellVariable hasOld() should set self.old
    101 Grids should take Lx, Ly, Lz arguments
    145 tests should be run with fipy.tests()
    177 remove ones and zeros from numerix.py
    178 Default time steps should be infinite
    291 term multiplication changes result
    296 FAQ gives bad guidance for anisotropic diffusion
    297 Use physical velocity in the manual/FAQ
    298 mesh manipulation of periodic meshes leads to errors
    299 Give helpfull error on - or / of meshes
    301 wrong cell to cell normal in periodic meshes
    302 gnuplot1d gives error on plot of facevariable
    309 pypi is failing
    312 Fresh FiPy gives ""ImportError: No viewers found"""
    314 Absence of enthought.tvtk causes test failures
    319 mesh in FiPy name space
    324 --pysparse configuration should never attempt MPI imports
    327 factoryMeshes.py not up to date with respect to keyword arguments
    331 changed constraints don't propagate
    332 anisotropic diffusion and constraints don't mix
    333 `--Trilinos --no-pysparse` uses PySparse?!?
    336 Profile and merge reconstrain branch
    339 close out reconstrain branch
    341 Fix fipy.terms._BinaryTerm test failure in parallel
    343 diffusionTerm(var=var1).solver(var=var0) should fail sensibly
    346 TeX is wrong in examples.phase.quaternary
    348 Include Benny's improved interpolation patch
    354 GmshExport is not tested and does not work
    355 Introduce mesh.x as shorthand for mesh.cellCenters[0] etc
    356 GmshImport should support all element types
    357 GmshImport should read element colors
    363 Reduce the run times for chemotaxis tests
    366 tests take *too* long!!!
    369 Make DiffusionTermNoCorrection the default
    370 Epetra Norm2 failure in parallel
    373 remove deprecated `steps=` from Solver
    376 remove deprecated `diffusionTerm=` argument to ConvectionTerm
    377 remove deprecated `NthOrderDiffusionTerm`
    380 remove deprecated Variable.transpose()
    381 remove deprecated viewers.make()
    382 get running in Py3k
    384 gmsh importer and gmsh tests don't clean up after themselves
    385 `diffusionTerm._test()` requires PySparse
    390 Improve test reporting to avoid inconsequential buildbot failures
    391 efficiency_test chokes on liquidVapor2D.py
    393 two `--scipy` failures
    395 `--pysparse --inline` failures
    417 Memory consumption growth with repeated meshing, especially with Gmsh
    418 Viewers not working when plotting meshes with zero cells in parallel
    419 examples/cahnHilliard/mesh2D.py broken with --trilinos
    420 Epetra.PyComm() broken on Debian
    421 cellVariable.min() broken in parallel
    426 Add in parallel buildbot testing on more than 2 processors
    427 Slow PyAMG solutions
    434 Gmsh I/O
    438 changes to gmshImport.py caused --inline problems
    439 gmshImport tests fail on Windows due to shared file
    441 Explicit convetion terms should fail when the equation has no TransientTerm (dt=None)
    445 getFaceCenters() should return a FaceVariable
    446 constraining values with ImplictSourceTerm not documented?
    448 Gmsh2D does not respect background mesh
    452 Gmsh background mesh doesn't work in parallel
    453 faceValue as FaceCenters gives inline failures
    454 Py3k and Windows test failures

.. warning::

   :term:`FiPy` 3 brought unavoidable syntax changes from :term:`FiPy` 2.
   Please see :mod:`examples.updating.update2_0to3_0` for guidance on the
   changes that you will need to make to your :term:`FiPy` 2.x scripts.

-------------------------
Download and Installation
-------------------------

Please refer to :ref:`INSTALLATION` for details on download and
installation. :term:`FiPy` can be redistributed and/or modified
freely, provided that any derivative works bear some notice that they
are derived from it, and any modified versions bear some notice that
they have been modified.

-------
Support
-------

You can communicate with the :term:`FiPy` developers and with other
users via our `mailing list`_ and we welcome you to use the `tracking
system`_ for bugs, support requests, feature requests and patch
submissions <http://matforge.org/fipy/report>. We welcome
collaborative efforts on this project.

:term:`FiPy` is a member of MatForge_, a project of the `Materials
Digital Library Pathway`_. This National Science Foundation funded
service provides management of our public source code repository, our
bug tracking system, and a "wiki" space for public contributions of
code snippets, discussions, and tutorials.

.. toctree::

   documentation/MAIL

------------------------
Conventions and Notation
------------------------

:term:`FiPy` is driven by :term:`Python` script files than you can view or modify in any
text editor.  :term:`FiPy` sessions are invoked from a command-line shell, such
as :command:`tcsh` or :command:`bash`.

Throughout, text to be typed at the keyboard will appear ``like this``.
Commands to be issued from an interactive shell will appear::

    $ like this

where you would enter the text ("``like this``") following the shell prompt,
denoted by "``$``".

Text blocks of the form::

    >>> a = 3 * 4
    >>> a
    12
    >>> if a == 12:
    ...     print "a is twelve"
    ...
    a is twelve

are intended to indicate an interactive session in the :term:`Python` interpreter.
We will refer to these as "interactive sessions" or as "doctest blocks".
The text "``>>>``" at the beginning of a line denotes the *primary prompt*,
calling for input of a :term:`Python` command.  The text "``...``" denotes the
*secondary prompt*, which calls for input that continues from the line
above, when required by :term:`Python` syntax.  All remaining lines, which begin
at the left margin, denote output from the :term:`Python` interpreter.  In all
cases, the prompt is supplied by the :term:`Python` interpreter and should not be
typed by you.

.. warning::

   :term:`Python` is sensitive to indentation and care should be taken to enter
   text exactly as it appears in the examples.

When references are made to file system paths, it is assumed that the
current working directory is the :term:`FiPy` distribution directory, refered to
as the "base directory", such that::

    examples/diffusion/steadyState/mesh1D.py

will correspond to, *e.g.*::

    /some/where/FiPy-X.Y/examples/diffusion/steadyState/mesh1D.py

Paths will always be rendered using POSIX conventions (path elements
separated by "``/``").  Any references of the form::

    examples.diffusion.steadyState.mesh1D

are in the :term:`Python` module notation and correspond to the equivalent POSIX
path given above.

We may at times use a 

.. note::

   to indicate something that may be of interest

or a

.. warning::

   to indicate something that could cause serious problems.

.. _MML:                  http://www.nist.gov/mml/
.. _CTCMS:                http://www.ctcms.nist.gov/
.. _MSED:                 http://www.nist.gov/mml/msed/
.. _NIST:                 http://www.nist.gov/
.. _compressed archive:   http://www.ctcms.nist.gov/fipy/download/FiPy-1.1.tar.gz
.. _tracking system:      http://matforge.org/fipy/report
.. _mailing list:         http://www.ctcms.nist.gov/fipy/documentation/MAIL.html
.. _Sourceforge:          http://www.sourceforge.net/projects/fipy
.. _Materials Digital Library Pathway: http://matdl.org
.. _MatForge:             http://matforge.org/

