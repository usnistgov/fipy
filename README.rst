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

| |TravisCI|_ |AppVeyor|_
| |GitHub|_ |gitter|_ |PyPi|_  |Codacy|_ |Depsy|_ |OpenHub|_ |CondaForge|_ |Binder|_

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

The significant changes since version 3.0 are:

- Level sets are now handled by :ref:`LSMLIBDOC` or :ref:`SCIKITFMM`
  solver libraries. These libraries are orders of magnitude faster than the
  original, :term:`Python`-only prototype.
- The :term:`Matplotlib` :func:`streamplot()` function can be used to display
  vector fields.
- Version control was switched to the Git_ distributed version control
  system. This system should make it much easier for :term:`FiPy` users to
  participate in development.

Tickets fixed in this release::

    62  "Move 'ImplicitDiffusionTerm().solve(var) == 0' ""failure"" from examples.phase.simple to examples.diffusion.mesh1D?"
    118 subscriber()._markStale() AttributeError
    138 `numerix.dot` doesn't support tensors
    143 "Trying to ""solve"" an integer `CellVariable` should raise an error"
    195 broken arithmetic face to cell distance calculations
    197 ~binOp doesn't work on branches/version-2_0
    305 add rhie chow correction term in stokes cavity example
    321 Windows interactive plotting mostly broken
    324 --pysparse configuration should never attempt MPI imports
    341 Fix fipy.terms._BinaryTerm test failure in parallel
    365 Rename GridXD
    368 Error adding meshes
    370 Epetra Norm2 failure in parallel
    383 move FiPy to distributed version control
    385 `diffusionTerm._test()` requires PySparse
    391 efficiency_test chokes on liquidVapor2D.py
    432 LSMLIB refactor
    441 Explicit convetion terms should fail when the equation has no TransientTerm (dt=None)
    445 getFaceCenters() should return a FaceVariable
    448 Gmsh2D does not respect background mesh
    452 Gmsh background mesh doesn't work in parallel
    453 faceValue as FaceCenters gives inline failures
    454 Assorted errors
    456 Web page links seem to be broken
    457 Make the citation links go to the DOI links
    460 Clean up interaction between dependencies and installation process
    461 SvnToGit clean up
    462 Fix for test failures on loki
    465 sign issues for equation with transient, convection and implicit terms
    466 "multiplying equation by ""x"" changes the solution"
    469 text in source:trunk/examples/convection/source.py is out of date
    470 Include mailing list activity frame on front page
    473 Gmsh importer can't read mesh elements with no tags
    475 getVersion() fails on Py3k
    477 Update Ohloh to point at git repo
    480 link to mailing list is wrong
    481 constrain should return a handle to the constraint for later deletion
    484 NIST CSS changed
    486 Using `Popen('gmsh ...', shell=True)` rather than `shell=False` security danger
    490 Parallel bug in non-uniform grids and conflicting mesh class and factory function names
    491 Rename communicator instances
    492 unOps can't be pickled
    493 Change documentation to promote use of stackoverflow
    494 Viewers don't inline well in IPython notebook
    496 FIPY_DISPLAY_MATRIX is broken
    497 examples/phase/binary.py has problems
    513 convection problem with cylindrical grid
    539 Bug with numpy 1.7.0
    557 NumPy 1.7.0 doesn't have _formatInteger
    564 VanLeerConvectionTerm MinMod slope limiter is broken
    638 numpy 1.7.1 test failures with physicalField.py
    639 Neumann boundary conditions not clearly documented
    641 Add support for Matplotlib streamplot
    648 Peclet inequalities have the wrong sign
    650 CylindricalNonUniformGrid2D doesn't make a FaceVariable for exteriorFaces
    652 Documentation change for Ubuntu install
    653 enable google analytics
    654 Switch to sphinxcontrib-bibtex
    655 Home page needs out-of-NIST redirects

.. warning::

   :term:`FiPy` 3 brought unavoidable syntax changes from :term:`FiPy` 2.
   Please see :mod:`examples.updating.update2_0to3_0` for guidance on the
   changes that you will need to make to your :term:`FiPy` 2.x scripts.

.. _Git: http://git-scm.com/

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
users via our `mailing list`_ and we welcome you to use the `issue
tracker`_ for bugs, support requests, feature requests and patch
submissions <https://github.com/usnistgov/fipy/issues>. We also monitor
StackOverflow_ for questions tagged with "fipy". We welcome
collaborative efforts on this project.

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

.. _MML:           http://www.nist.gov/mml/
.. _CTCMS:         http://www.ctcms.nist.gov/
.. _MSED:          http://www.nist.gov/mml/msed/
.. _NIST:          http://www.nist.gov/
.. _issue tracker: https://github.com/usnistgov/fipy/issues
.. _mailing list:  http://www.ctcms.nist.gov/fipy/documentation/MAIL.html
.. _StackOverflow: http://stackoverflow.com/questions/tagged/fipy

.. |GitHub|        image:: https://img.shields.io/github/contributors/usnistgov/fipy.svg
.. _Github:        https://github.com/usnistgov/fipy
.. |gitter|        image:: https://badges.gitter.im/usnistgov/fipy.svg
.. _gitter:        https://gitter.im/usnistgov/fipy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=body_badge
.. |TravisCI|      image:: https://img.shields.io/travis/usnistgov/fipy/develop.svg?label=Linux
.. _TravisCI:      https://travis-ci.org/usnistgov/fipy
.. |AppVeyor|      image:: https://ci.appveyor.com/api/projects/status/github/usnistgov/fipy?branch=develop&svg=true&failingText=Windows%20-%20failing&passingText=Windows%20-%20passing
.. _AppVeyor:      https://ci.appveyor.com/project/guyer/fipy
.. |OpenHub|       image:: https://www.openhub.net/p/fipy/widgets/project_thin_badge.gif
.. _OpenHub:       https://www.openhub.net/p/fipy
.. |PyPi|          image:: https://img.shields.io/pypi/v/fipy.svg
.. _PyPi:          https://pypi.python.org/pypi/FiPy
.. |CondaForge|    image:: https://anaconda.org/guyer/fipy/badges/downloads.svg
.. _CondaForge:    https://anaconda.org/guyer/fipy
.. |Depsy|         image:: http://depsy.org/api/package/pypi/FiPy/badge.svg
.. _Depsy:         http://depsy.org/package/python/FiPy
.. |Codacy|         image:: https://api.codacy.com/project/badge/Grade/d02921bb54b14e88a1e2e1f5520133f4
.. _Codacy:         https://www.codacy.com/app/tkphd/fipy?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=usnistgov/fipy&amp;utm_campaign=Badge_Grade
.. |Binder|        image:: https://mybinder.org/badge.svg
.. _Binder:        https://mybinder.org/v2/gh/usnistgov/fipy/develop?filepath=examples%2Findex.ipynb
