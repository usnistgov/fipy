.. Generate incremental updates to this file with
   $ python setup.py changelog <OPTIONS>

.. _CHANGELOG:

==========
Change Log
==========

-----------------
Version |release|
-----------------

This maintenance release:

- Addresses compatibility with recent releases of
  Python 3.12, NumPy 2.0, SciPy 1.14, and PETSc 3.20.
- Adds `conda-lock <https://github.com/conda/conda-lock>`_ environment
  lock files with specified compatible versions of FiPy prerequisites.
- Fixes numerous documentation errors.

.. attention::

   SciPy 1.13.0 generates one test suite error for
   ``fipy.matrices.scipyMatrix._ScipyMatrix.CSR``.  Either ignore the test
   failure or upgrade to SciPy >= 1.13.1

.. attention::

   PETSc 3.21 crashes our test suite when running in parallel (`#1054
   <https://github.com/usnistgov/fipy/issues/1054>`_).  PETSc <= 3.20 is
   recommended, although `petsc 3.20.2_*_102 is broken on macOS
   <https://github.com/conda-forge/petsc-feedstock/issues/180>`_.

Pulls
-----

- Introduce Timer context manager
  (`#995 <https://github.com/usnistgov/fipy/pull/995>`_)
- switch nix recipe to flake
  (`#992 <https://github.com/usnistgov/fipy/pull/992>`_)
- Tweak documentation
  (`#991 <https://github.com/usnistgov/fipy/pull/991>`_)
- Log much more information about FiPy environment
  (`#990 <https://github.com/usnistgov/fipy/pull/990>`_)
- Fix inclusion of `environments/README.rst`
  (`#988 <https://github.com/usnistgov/fipy/pull/988>`_)
- Environment pinning
  (`#985 <https://github.com/usnistgov/fipy/pull/985>`_)

Fixes
-----

- `#1049 <https://github.com/usnistgov/fipy/issues/1049>`_:
  Numpy 2.0.0 breaks things
- `#1010 <https://github.com/usnistgov/fipy/issues/1010>`_:
  `examples.diffusion.mesh1D` No-flux - steady-state doesn't always give
  zero
- `#1000 <https://github.com/usnistgov/fipy/issues/1000>`_:
  `examples.diffusion.mesh1D` constrains a gradient but calls it a flux
- `#997 <https://github.com/usnistgov/fipy/issues/997>`_:
  `future.standard_library` breaking python 3.12 compatibility
- `#967 <https://github.com/usnistgov/fipy/issues/967>`_:
  Sign error in Robin condition
- `#963 <https://github.com/usnistgov/fipy/issues/963>`_:
  PETSc 3.20.0 broke the world
- `#961 <https://github.com/usnistgov/fipy/issues/961>`_:
  Representation of index variables is broken
- `#952 <https://github.com/usnistgov/fipy/issues/952>`_:
  Uncaught Exception from the no-flux steady-state diffusion example
- `#944 <https://github.com/usnistgov/fipy/issues/944>`_:
  Having problem with Viewer
- `#865 <https://github.com/usnistgov/fipy/issues/865>`_:
  Sphinx search is broken on website
- `#673 <https://github.com/usnistgov/fipy/issues/673>`_:
  Deprecations don't properly format properties
- `#512 <https://github.com/usnistgov/fipy/issues/512>`_:
  Default coefficient of `ImplicitSourceTerm` is 0

--------------------------
Version 3.4.4 - 2023-06-27
--------------------------

This maintenance release adds :ref:`LOGGING` and resolves compatibility issues
with recent builds of :term:`PETSc` and :term:`NumPy`.

Pulls
-----

- Fix numpy 1.25 issues
  (`#930 <https://github.com/usnistgov/fipy/pull/930>`_)
- Get CI working again
  (`#925 <https://github.com/usnistgov/fipy/pull/925>`_)
- Discourage StackOverflow
  (`#876 <https://github.com/usnistgov/fipy/pull/876>`_)
- Add Logging
  (`#875 <https://github.com/usnistgov/fipy/pull/875>`_)
- Add tests for the Nix build
  (`#791 <https://github.com/usnistgov/fipy/pull/791>`_)

Fixes
-----

- `#896 <https://github.com/usnistgov/fipy/issues/896>`_:
  Poor garbage collection with petsc4py 3.18.3 (was "Memory leak in
  `term.justErrorVector()`", but this isn't strictly a leak)

--------------------------
Version 3.4.3 - 2022-06-15
--------------------------

This maintenance release adds a new example contributed by
`@Jon83Carvalho <https://github.com/Jon83Carvalho>`_,
clarifies many points in the documentation,
migrates all :ref:`CONTINUOUSINTEGRATION` to
`Azure <https://dev.azure.com>`_,
updates to using
`wheels <https://packaging.python.org/en/latest/specifications/binary-distribution-format/>`_
for distribution,
and substantially refactors matrices to work more consistently across
solvers.

Pulls
-----

- Update CI documentation to refer only to Azure
  (`#863 <https://github.com/usnistgov/fipy/pull/863>`_)
- Refine azure runs
  (`#851 <https://github.com/usnistgov/fipy/pull/851>`_)
- Debug CIs
  (`#848 <https://github.com/usnistgov/fipy/pull/848>`_)
- Collect contact information on single page
  (`#847 <https://github.com/usnistgov/fipy/pull/847>`_)
- Set up CI with Azure Pipelines
  (`#822 <https://github.com/usnistgov/fipy/pull/822>`_)
- Replace deprecated numpy types
  (`#798 <https://github.com/usnistgov/fipy/pull/798>`_)
- Move trilinos tests to Py3k
  (`#797 <https://github.com/usnistgov/fipy/pull/797>`_)
- Fix Python 2.7 conda environment
  (`#795 <https://github.com/usnistgov/fipy/pull/795>`_)
- fix: stop divide by zero warning in LU solvers
  (`#790 <https://github.com/usnistgov/fipy/pull/790>`_)
- Introduce `SharedTemporaryFile` (bis)
  (`#769 <https://github.com/usnistgov/fipy/pull/769>`_)
- Raise `ImportError` before trying to unpack solvers
  (`#768 <https://github.com/usnistgov/fipy/pull/768>`_)
- Disable TVTK tests if its prerequisites aren't met
  (`#764 <https://github.com/usnistgov/fipy/pull/764>`_)
- Tabulate versions of FiPy dependencies when tests are run
  (`#763 <https://github.com/usnistgov/fipy/pull/763>`_)
- Debug CI failures
  (`#749 <https://github.com/usnistgov/fipy/pull/749>`_)
- Stokes Cavity - non-Newtonian
  (`#748 <https://github.com/usnistgov/fipy/pull/748>`_)
  Thanks to `@Jon83Carvalho <https://github.com/Jon83Carvalho>`_.
- Refactor matrices
  (`#721 <https://github.com/usnistgov/fipy/pull/721>`_)

Fixes
-----

- `#862 <https://github.com/usnistgov/fipy/issues/862>`_:
  Could not load the Qt platform plugin "`xcb`"
- `#858 <https://github.com/usnistgov/fipy/issues/858>`_:
  CI issues
- `#856 <https://github.com/usnistgov/fipy/issues/856>`_:
  `FaceVariable` does not accumulate properly in parallel
- `#850 <https://github.com/usnistgov/fipy/issues/850>`_:
  Switch to wheels
- `#849 <https://github.com/usnistgov/fipy/issues/849>`_:
  `linux-py27-pysparse` fails
- `#841 <https://github.com/usnistgov/fipy/issues/841>`_:
  `Matplotlib2DViewer` should accept color map as string
- `#836 <https://github.com/usnistgov/fipy/issues/836>`_:
  Document that coupled and high-order diffusion terms are
  incompatible
- `#833 <https://github.com/usnistgov/fipy/issues/833>`_:
  `fipy.tools.dump` undocumented that it always gzips
- `#828 <https://github.com/usnistgov/fipy/issues/828>`_:
  `colorbar=True` no longer works Stokes flow example
- `#826 <https://github.com/usnistgov/fipy/issues/826>`_:
  Gmsh load issue
- `#818 <https://github.com/usnistgov/fipy/issues/818>`_:
  Document that `GridND` meshes are always Cartesian
- `#811 <https://github.com/usnistgov/fipy/issues/811>`_:
  In python 3.9 __repr__ throws an exception with abs
- `#801 <https://github.com/usnistgov/fipy/issues/801>`_:
  CircleCI test-36-trilinos-serial extremely slow
- `#800 <https://github.com/usnistgov/fipy/issues/800>`_:
  CircleCI conda2_env is really slow and ends up installing FiPy 3.3
- `#796 <https://github.com/usnistgov/fipy/issues/796>`_:
  `examples.phase.polyxtal` freezes on CircleCI with Py3k and scipy
  solvers
- `#792 <https://github.com/usnistgov/fipy/issues/792>`_:
  `circleQuad` example fails with Gmsh > 4.4
- `#781 <https://github.com/usnistgov/fipy/issues/781>`_:
  `MatplolibViewer.axes` property is not documented
- `#778 <https://github.com/usnistgov/fipy/issues/778>`_:
  Binder failed build
- `#762 <https://github.com/usnistgov/fipy/issues/762>`_:
  Equations on Website don't show right
- `#742 <https://github.com/usnistgov/fipy/issues/742>`_:
  No documentation for `Variable.mag`
- `#735 <https://github.com/usnistgov/fipy/issues/735>`_:
  `pip install fipy` fails
- `#734 <https://github.com/usnistgov/fipy/issues/734>`_:
  Document the residual
- `#688 <https://github.com/usnistgov/fipy/issues/688>`_:
  try-except not needed for circle Viewer
- `#676 <https://github.com/usnistgov/fipy/issues/676>`_:
  Default no-flux condition is not explicitly stated
- `#609 <https://github.com/usnistgov/fipy/issues/609>`_:
  Parallelizing of Gmsh meshes not clearly documented
- `#400 <https://github.com/usnistgov/fipy/issues/400>`_:
  Fix `FaceVariable.globalValue` method

----------------------------
Version 3.4.2.1 - 2020-08-01
----------------------------

This release fixes assorted viewer issues, fixes a problem with convection
boundary conditions, and introduces spherical meshes.

.. attention::

   There are
   `known <https://travis-ci.com/github/usnistgov/fipy/builds/177879719>`_
   `failures <https://app.circleci.com/pipelines/github/usnistgov/fipy/248/workflows/4babcd98-aafc-4931-a353-64bbb3c93cb6>`_
   with the VTK viewers (bitrot has started to set
   in since the `demise of Python 2.7`_).  There's also a new parallel
   failure in `NonUniformGrid1D` that we need to figure out.

.. _demise of Python 2.7: https://www.python.org/dev/peps/pep-0373/#update

Pulls
-----

- Move mailing list
  (`#747 <https://github.com/usnistgov/fipy/pull/747>`_)
- `Spherical1D` (`Uniform` and `NonUniform`) meshes
  (`#732 <https://github.com/usnistgov/fipy/pull/732>`_)
  Thanks to `@klkuhlm <https://github.com/klkuhlm>`_.
- fix Neumann BCs using constraints with convection terms
  (`#719 <https://github.com/usnistgov/fipy/pull/719>`_)
  Thanks to `@atismer <https://github.com/atismer>`_.
- Add vertex index inversions
  (`#716 <https://github.com/usnistgov/fipy/pull/716>`_)

Fixes
-----

- `#726 <https://github.com/usnistgov/fipy/issues/726>`_:
  `MayaviClient` not compatible with Python 3
- `#663 <https://github.com/usnistgov/fipy/issues/663>`_:
  `datamin`/`datamax` argument ignored by viewer
- `#662 <https://github.com/usnistgov/fipy/issues/662>`_:
  Issues Scaling `Colorbar` with `Datamin` and `Datamax` `Args`

--------------------------
Version 3.4.1 - 2020-02-14
--------------------------

This release is primarily for compatibility with :mod:`numpy` 1.18.

Pulls
-----

- Fix documentation
  (`#711 <https://github.com/usnistgov/fipy/pull/711>`_)
- build(nix): fix broken plm_rsh_agent error
  (`#710 <https://github.com/usnistgov/fipy/pull/710>`_)
- CIs error on deprecation warning
  (`#708 <https://github.com/usnistgov/fipy/pull/708>`_)

Fixes
-----

- `#703 <https://github.com/usnistgov/fipy/issues/703>`_:
  FORTRAN array ordering is deprecated

------------------------
Version 3.4 - 2020-02-06
------------------------

This release adds support for the :term:`PETSc` solvers for
:ref:`PARALLEL`.

Pulls
-----

- Add support for PETSc solvers
  (`#701 <https://github.com/usnistgov/fipy/pull/701>`_)
- Assorted fixes while supporting PETSc
  (`#700 <https://github.com/usnistgov/fipy/pull/700>`_)
  - Fix print statements for Py3k
  - Resolve Gmsh issues
  - Dump only on processor 0
  - Only write `timetests` on processor 0
  - Fix conda-forge link
  - Upload PDF
  - Document `print` option of `FIPY_DISPLAY_MATRIX`
  - Use legacy numpy formatting when testing individual modules
  - Switch to matplotlib's built-in symlog scaling
  - Clean up tests
- Assorted fixes for benchmark 8
  (`#699 <https://github.com/usnistgov/fipy/pull/699>`_)
  - Stipulate `--force` option for `conda remove fipy`
  - Update Miniconda installation url
  - Replace `_CellVolumeAverageVariable` class with `Variable` expression
  - Fix output for bad call stack
- Make CircleCI build docs on Py3k
  (`#698 <https://github.com/usnistgov/fipy/pull/698>`_)
- Fix link to Nick Croft's thesis
  (`#681 <https://github.com/usnistgov/fipy/pull/681>`_)
- Fix NIST header footer
  (`#680 <https://github.com/usnistgov/fipy/pull/680>`_)
- Use Nixpkgs version of FiPy expression
  (`#661 <https://github.com/usnistgov/fipy/pull/661>`_)
- Update the Nix recipe
  (`#658 <https://github.com/usnistgov/fipy/pull/658>`_)

Fixes
-----

- `#692 <https://github.com/usnistgov/fipy/issues/692>`_:
  Can't copy example scripts with the command line
- `#669 <https://github.com/usnistgov/fipy/issues/669>`_:
  input() deadlock on parallel runs
- `#643 <https://github.com/usnistgov/fipy/issues/643>`_:
  Automate release process

------------------------
Version 3.3 - 2019-06-28
------------------------

This release brings support for Python 2 and Python 3 from the same source,
without any translation.  Thanks to `@pya <https://github.com/pya>`_ and
`@woodscn <https://github.com/woodscn>`_ for getting things started.

Pulls
-----

- Automate spell check
  (`#657 <https://github.com/usnistgov/fipy/pull/657>`_)
- Fix gmsh on windows
  (`#648 <https://github.com/usnistgov/fipy/pull/648>`_)
- Fix sphinx documentation
  (`#647 <https://github.com/usnistgov/fipy/pull/647>`_)
- Migrate to Py3k
  (`#645 <https://github.com/usnistgov/fipy/pull/645>`_)
- `gmshMesh.py` compatibility with Gmsh > 3.0.6
  (`#644 <https://github.com/usnistgov/fipy/pull/644>`_)
  Thanks to `@xfong <https://github.com/xfong>`_.

Fixes
-----

- `#655 <https://github.com/usnistgov/fipy/issues/655>`_:
  When Python 2 and 3 are installed, Mayavi wont work.
  Thanks to `@Hendrik410 <https://github.com/Hendrik410>`_.
- `#646 <https://github.com/usnistgov/fipy/issues/646>`_:
  Deprecate develop branch
- `#643 <https://github.com/usnistgov/fipy/issues/643>`_:
  Automate release process
- `#601 <https://github.com/usnistgov/fipy/issues/601>`_:
  :file:`contents.rst` and :file:`manual.rst` are a recursive mess
- `#597 <https://github.com/usnistgov/fipy/issues/597>`_:
  Use GitHub link for the compressed archive in documentation
- `#557 <https://github.com/usnistgov/fipy/issues/557>`_:
  `faceGradAverage` is stupid
- `#552 <https://github.com/usnistgov/fipy/issues/552>`_:
  documentation integration
- `#458 <https://github.com/usnistgov/fipy/issues/458>`_:
  Documentation wrong for precedence of `Lx` and `dx` for
  `NonUniformGrids`
- `#457 <https://github.com/usnistgov/fipy/issues/457>`_:
  Special methods are not included in Sphinx documentation
- `#432 <https://github.com/usnistgov/fipy/issues/432>`_:
  Python 3 issues
- `#340 <https://github.com/usnistgov/fipy/issues/340>`_:
  Don't upload packages to PyPI, just add the master url

------------------------
Version 3.2 - 2019-04-22
------------------------

This is predominantly a `DevOps`_ release.  The focus has been on making
FiPy easier to install with :term:`conda`.  It's also possible to install a
minimal set of prerequisites with :term:`pip`.  Further, :term:`FiPy` is
automatically tested on all major platforms using cloud-based
:ref:`CONTINUOUSINTEGRATION` (:term:`linux` with :term:`CircleCI`,
:term:`macOS` with :term:`TravisCI`, and :term:`Windows` with
:term:`AppVeyor`).

Pulls
-----

- Make badges work in GitHub and pdf
  (`#636 <https://github.com/usnistgov/fipy/pull/636>`_)
- Fix Robin errors
  (`#615 <https://github.com/usnistgov/fipy/pull/615>`_)
- Issue555 inclusive license
  (`#613 <https://github.com/usnistgov/fipy/pull/613>`_)
- Update CIs
  (`#607 <https://github.com/usnistgov/fipy/pull/607>`_)
- Add CHANGELOG and tool to generate from issues and pull requests
  (`#600 <https://github.com/usnistgov/fipy/pull/600>`_)
- Explain where to get examples
  (`#596 <https://github.com/usnistgov/fipy/pull/596>`_)
- spelling corrections using en_US dictionary
  (`#594 <https://github.com/usnistgov/fipy/pull/594>`_)
- Remove `SmoothedAggregationSolver`
  (`#593 <https://github.com/usnistgov/fipy/pull/593>`_)
- Nix recipe for FiPy
  (`#585 <https://github.com/usnistgov/fipy/pull/585>`_)
- Point PyPI to github master tarball
  (`#582 <https://github.com/usnistgov/fipy/pull/582>`_)
- Revise Navier-Stokes expression in the viscous limit
  (`#580 <https://github.com/usnistgov/fipy/pull/580>`_)
- Update `stokesCavity.py`
  (`#579 <https://github.com/usnistgov/fipy/pull/579>`_)
  Thanks to `@Rowin <https://github.com/Rowin>`_.
- Add `--inline` to TravisCI tests
  (`#578 <https://github.com/usnistgov/fipy/pull/578>`_)
- Add support for binder
  (`#577 <https://github.com/usnistgov/fipy/pull/577>`_)
- Fix `epetra vector not numarray`
  (`#574 <https://github.com/usnistgov/fipy/pull/574>`_)
- add Codacy badge
  (`#572 <https://github.com/usnistgov/fipy/pull/572>`_)
- Fix output when PyTrilinos or PyTrilinos version is unavailable
  (`#570 <https://github.com/usnistgov/fipy/pull/570>`_)
  Thanks to `@shwina <https://github.com/shwina>`_.
- Fix check for PyTrilinos
  (`#569 <https://github.com/usnistgov/fipy/pull/569>`_)
  Thanks to `@shwina <https://github.com/shwina>`_.
- Adding support for GPU solvers via pyamgx
  (`#567 <https://github.com/usnistgov/fipy/pull/567>`_)
  Thanks to `@shwina <https://github.com/shwina>`_.
- revise dedication to the public domain
  (`#556 <https://github.com/usnistgov/fipy/pull/556>`_)
- Fix tests that don't work in parallel
  (`#550 <https://github.com/usnistgov/fipy/pull/550>`_)
- add badges to index and readme
  (`#546 <https://github.com/usnistgov/fipy/pull/546>`_)
- Ensure vector is `dtype` float before matrix multiply
  (`#544 <https://github.com/usnistgov/fipy/pull/544>`_)
- Revert "Issue534 physical field mishandles compound units"
  (`#536 <https://github.com/usnistgov/fipy/pull/536>`_)
- Document boundary conditions
  (`#532 <https://github.com/usnistgov/fipy/pull/532>`_)
- Deadlocks and races
  (`#524 <https://github.com/usnistgov/fipy/pull/524>`_)
- Make max/min global
  (`#520 <https://github.com/usnistgov/fipy/pull/520>`_)
- Add a Gitter chat badge to :file:`README.rst`
  (`#516 <https://github.com/usnistgov/fipy/pull/516>`_)
  Thanks to `@gitter-badger <https://github.com/gitter-badger>`_.
- Add TravisCI build recipe
  (`#489 <https://github.com/usnistgov/fipy/pull/489>`_)

Fixes
-----

- `#631 <https://github.com/usnistgov/fipy/issues/631>`_:
  Clean up :file:`INSTALLATION.rst`
- `#628 <https://github.com/usnistgov/fipy/issues/628>`_:
  Problems with the viewer
- `#627 <https://github.com/usnistgov/fipy/issues/627>`_:
  Document OMP_NUM_THREADS
- `#625 <https://github.com/usnistgov/fipy/issues/625>`_:
  `setup.py` should not import fipy
- `#623 <https://github.com/usnistgov/fipy/issues/623>`_:
  Start using `versioneer`
- `#621 <https://github.com/usnistgov/fipy/issues/621>`_:
  Plot `FaceVariable` with matplotlib
- `#617 <https://github.com/usnistgov/fipy/issues/617>`_:
  Pick 1st Value and last Value of 1D `CellVariable` while running in
  parallel
- `#611 <https://github.com/usnistgov/fipy/issues/611>`_:
  The coefficient cannot be a `FaceVariable` ??
- `#610 <https://github.com/usnistgov/fipy/issues/610>`_:
  Anisotropy example: Contour plot displaying in legend of figure !?
- `#608 <https://github.com/usnistgov/fipy/issues/608>`_:
  `var.mesh`: `Property` object not callable...?
- `#603 <https://github.com/usnistgov/fipy/issues/603>`_:
  Can't run basic test or examples
- `#602 <https://github.com/usnistgov/fipy/issues/602>`_:
  Revise build and release documentation
- `#592 <https://github.com/usnistgov/fipy/issues/592>`_:
  is :file:`resources.rst` useful?
- `#590 <https://github.com/usnistgov/fipy/issues/590>`_:
  No module named `pyAMGSolver`
- `#584 <https://github.com/usnistgov/fipy/issues/584>`_:
  Viewers don't animate in jupyter notebook
- `#566 <https://github.com/usnistgov/fipy/issues/566>`_:
  Support for GPU solvers using pyamgx
- `#565 <https://github.com/usnistgov/fipy/issues/565>`_:
  pip install does not work on empty env
- `#564 <https://github.com/usnistgov/fipy/issues/564>`_:
  Get green boxes across the board
- `#561 <https://github.com/usnistgov/fipy/issues/561>`_:
  Cannot cast array data from `dtype('int64')` to `dtype('int32')`
  according to the rule `safe`
- `#555 <https://github.com/usnistgov/fipy/issues/555>`_:
  inclusive license
- `#551 <https://github.com/usnistgov/fipy/issues/551>`_:
  Sphinx spews many warnings:
- `#545 <https://github.com/usnistgov/fipy/issues/545>`_:
  Many Py3k failures
- `#543 <https://github.com/usnistgov/fipy/issues/543>`_:
  Epetra Vector can't be integer
- `#539 <https://github.com/usnistgov/fipy/issues/539>`_:
  `examples/diffusion/explicit/mixedElement.py` is a mess
- `#538 <https://github.com/usnistgov/fipy/issues/538>`_:
  badges
- `#534 <https://github.com/usnistgov/fipy/issues/534>`_:
  `PhysicalField` mishandles compound units
- `#533 <https://github.com/usnistgov/fipy/issues/533>`_:
  pip or conda installation don't make clear where to get examples
- `#531 <https://github.com/usnistgov/fipy/issues/531>`_:
  `drop_tol` argument to `scipy.sparse.linalg.splu` is gone
- `#530 <https://github.com/usnistgov/fipy/issues/530>`_:
  conda installation instructions not explicit about python version
- `#528 <https://github.com/usnistgov/fipy/issues/528>`_:
  scipy 1.0.0 incompatibilities
- `#525 <https://github.com/usnistgov/fipy/issues/525>`_:
  conda `guyer/pysparse` doesn't run on osx
- `#513 <https://github.com/usnistgov/fipy/issues/513>`_:
  Stokes example gives wrong equation
- `#510 <https://github.com/usnistgov/fipy/issues/510>`_:
  Weave, Scipy and `--inline`
- `#509 <https://github.com/usnistgov/fipy/issues/509>`_:
  Unable to use conda for installing FiPy in Windows
- `#506 <https://github.com/usnistgov/fipy/issues/506>`_:
  Error using spatially varying anisotropic diffusion coefficient
- `#488 <https://github.com/usnistgov/fipy/issues/488>`_:
  Gmsh 2.11 breaks `GmshGrids`
- `#435 <https://github.com/usnistgov/fipy/issues/435>`_:
  `pip install pysparse` fails with
  "`fatal error: 'spmatrix.h' file not found`"
- `#434 <https://github.com/usnistgov/fipy/issues/434>`_:
  `pip install fipy` fails with 
  "`ImportError: No module named ez_setup`"

.. _DevOps:   https://en.wikipedia.org/wiki/DevOps

--------------------------
Version 3.1.3 - 2017-01-17
--------------------------

Fixes
-----

- `#502 <https://github.com/usnistgov/fipy/issues/502>`_:
  gmane is defunct

--------------------------
Version 3.1.2 - 2016-12-24
--------------------------

Pulls
-----

- remove `recvobj` from calls to `allgather`, require `sendobj`
  (`#492 <https://github.com/usnistgov/fipy/pull/492>`_)
- restore trailing whitespace to expected output of pysparse matrix
  tests
  (`#485 <https://github.com/usnistgov/fipy/pull/485>`_)
- Format version string for pep 440
  (`#483 <https://github.com/usnistgov/fipy/pull/483>`_)
- Provide some documentation for what `_faceToCellDistanceRatio` is
  and why it's scalar
  (`#481 <https://github.com/usnistgov/fipy/pull/481>`_)
- Strip all trailing white spaces and empty lines at EOF for `.py` and
  `.r`?
  (`#479 <https://github.com/usnistgov/fipy/pull/479>`_)
  Thanks to `@pya <https://github.com/pya>`_.
- `fipy/meshes/uniformGrid3D.py`: fix `_cellToCellIDs` and more
  `concatenate()` calls
  (`#478 <https://github.com/usnistgov/fipy/pull/478>`_)
  Thanks to `@pkgw <https://github.com/pkgw>`_.
- Remove incorrect `axis` argument to `concatenate`
  (`#477 <https://github.com/usnistgov/fipy/pull/477>`_)
- Updated to NumPy 1.10
  (`#472 <https://github.com/usnistgov/fipy/pull/472>`_)
  Thanks to `@pya <https://github.com/pya>`_.
- Some spelling corrections
  (`#471 <https://github.com/usnistgov/fipy/pull/471>`_)
  Thanks to `@pkgw <https://github.com/pkgw>`_.
- Sort entry points by package name before testing.
  (`#469 <https://github.com/usnistgov/fipy/pull/469>`_)
- Update import syntax in examples
  (`#466 <https://github.com/usnistgov/fipy/pull/466>`_)
- Update links to prerequisites
  (`#465 <https://github.com/usnistgov/fipy/pull/465>`_)
- Correct implementation of `examples.cahnHilliard.mesh2DCoupled`. Fixes
  ?
  (`#463 <https://github.com/usnistgov/fipy/pull/463>`_)
- Fix typeset analytical solution
  (`#460 <https://github.com/usnistgov/fipy/pull/460>`_)
- Clear `pdflatex` build errors by removing :term:`Python` from heading
  (`#459 <https://github.com/usnistgov/fipy/pull/459>`_)
- purge gist from viewers and optional module lists in `setup.py`
  (`#456 <https://github.com/usnistgov/fipy/pull/456>`_)
- Remove deprecated methods that duplicate NumPy ufuncs
  (`#454 <https://github.com/usnistgov/fipy/pull/454>`_)
- Remove deprecated Gmsh importers
  (`#452 <https://github.com/usnistgov/fipy/pull/452>`_)
- Remove deprecated getters and setters
  (`#450 <https://github.com/usnistgov/fipy/pull/450>`_)
- Update links for FiPy developers
  (`#448 <https://github.com/usnistgov/fipy/pull/448>`_)
- Render appropriately if in IPython notebook
  (`#447 <https://github.com/usnistgov/fipy/pull/447>`_)
- Plot contour in proper axes
  (`#446 <https://github.com/usnistgov/fipy/pull/446>`_)
- Robust Gmsh version checking with `distutils.version.StrictVersion`
  (`#442 <https://github.com/usnistgov/fipy/pull/442>`_)
- compare gmsh versions as tuples, not floats
  (`#441 <https://github.com/usnistgov/fipy/pull/441>`_)
- Corrected two tests
  (`#439 <https://github.com/usnistgov/fipy/pull/439>`_)
  Thanks to `@alfrenardi <https://github.com/alfrenardi>`_.
- Issue426 fix robin example typo
  (`#431 <https://github.com/usnistgov/fipy/pull/431>`_)
  Thanks to `@raybsmith <https://github.com/raybsmith>`_.
- Issue426 fix robin example analytical solution
  (`#429 <https://github.com/usnistgov/fipy/pull/429>`_)
  Thanks to `@raybsmith <https://github.com/raybsmith>`_.
- Force `MatplotlibViewer` to display
  (`#428 <https://github.com/usnistgov/fipy/pull/428>`_)
- Allow for 2 periodic axes in 3D
  (`#424 <https://github.com/usnistgov/fipy/pull/424>`_)
- Bug with Matplotlib 1.4.0 is fixed
  (`#419 <https://github.com/usnistgov/fipy/pull/419>`_)

Fixes
-----

- `#498 <https://github.com/usnistgov/fipy/issues/498>`_:
  nonlinear source term
- `#496 <https://github.com/usnistgov/fipy/issues/496>`_:
  `scipy.LinearBicgstabSolver` doesn't take arguments
- `#494 <https://github.com/usnistgov/fipy/issues/494>`_:
  Gmsh call errors
- `#493 <https://github.com/usnistgov/fipy/issues/493>`_:
  `Reviewable.io` has read-only access, can't leave comments
- `#491 <https://github.com/usnistgov/fipy/issues/491>`_:
  `globalValue` raises error from mpi4py
- `#484 <https://github.com/usnistgov/fipy/issues/484>`_:
  Pysparse tests fail
- `#482 <https://github.com/usnistgov/fipy/issues/482>`_:
  FiPy development version string not compliant with PEP 440
- `#476 <https://github.com/usnistgov/fipy/issues/476>`_:
  `setuptools` 18.4 breaks test suite
- `#475 <https://github.com/usnistgov/fipy/issues/475>`_:
  `Grid3D` broken by numpy 1.10
- `#470 <https://github.com/usnistgov/fipy/issues/470>`_:
  `Mesh3D` `cellToCellIDs` is broken
- `#467 <https://github.com/usnistgov/fipy/issues/467>`_:
  Out-of-sequence Viewer imports
- `#462 <https://github.com/usnistgov/fipy/issues/462>`_:
  GMSH version >= 2.10 incorrectly read by `gmshMesh.py`
- `#455 <https://github.com/usnistgov/fipy/issues/455>`_:
  `setup.py` gist warning
- `#445 <https://github.com/usnistgov/fipy/issues/445>`_:
  `DendriteViewer` puts contours over color bar
- `#443 <https://github.com/usnistgov/fipy/issues/443>`_:
  `MatplotlibViewer` still has problems in IPython notebook
- `#440 <https://github.com/usnistgov/fipy/issues/440>`_:
  Use github API to get nicely formatted list of issues
- `#438 <https://github.com/usnistgov/fipy/issues/438>`_:
  Failed tests on Mac OS X
- `#437 <https://github.com/usnistgov/fipy/issues/437>`_:
  Figure misleading in `examples.cahnHilliard.mesh2DCoupled`
- `#433 <https://github.com/usnistgov/fipy/issues/433>`_:
  Links to prerequisites are broken
- `#430 <https://github.com/usnistgov/fipy/issues/430>`_:
  Make develop the default branch on Github
- `#427 <https://github.com/usnistgov/fipy/issues/427>`_:
  `MatplotlibViewer` don't display
- `#425 <https://github.com/usnistgov/fipy/issues/425>`_:
  Links for Warren and Guyer are broken on the web page
- `#421 <https://github.com/usnistgov/fipy/issues/421>`_:
  The "limits" argument for `Matplotlib2DGridViewer` does not function
- `#416 <https://github.com/usnistgov/fipy/issues/416>`_:
  Updates to reflect move to Github

--------------------------
Version 3.1.1 - 2015-12-17
--------------------------

Fixes
-----

- `#415 <https://github.com/usnistgov/fipy/issues/415>`_:
  `MatplotlibGrid2DViewer` error with Matplotlib version 1.4.0
- `#414 <https://github.com/usnistgov/fipy/issues/414>`_:
  `PeriodicGrid3D` supports Only 1 axes of periodicity or all 3, not 2
- `#413 <https://github.com/usnistgov/fipy/issues/413>`_:
  Remind users of different types of conservation equations
- `#412 <https://github.com/usnistgov/fipy/issues/412>`_:
  Pickling Communicators is unnecessary for Grids
- `#408 <https://github.com/usnistgov/fipy/issues/408>`_:
  Implement `PeriodicGrid3D`
- `#407 <https://github.com/usnistgov/fipy/issues/407>`_:
  Strange deprecation loop in reshape()
- `#404 <https://github.com/usnistgov/fipy/issues/404>`_:
  package never gets uploaded to PyPI
- `#401 <https://github.com/usnistgov/fipy/issues/401>`_:
  Vector equations are broken when `sweep` is used instead of `solve`.
- `#295 <https://github.com/usnistgov/fipy/issues/295>`_:
  Gmsh version must be >= 2.0 errors on `zizou`

------------------------
Version 3.1 - 2013-09-30
------------------------

The significant changes since version 3.0 are:

- Level sets are now handled by :ref:`LSMLIBDOC` or :ref:`SCIKITFMM` 
  solver libraries. These libraries are orders of magnitude faster than the 
  original, :term:`Python`-only prototype.
- The :term:`Matplotlib` :func:`streamplot()` function can be used to display 
  vector fields.
- Version control was switched to the Git_ distributed version control 
  system. This system should make it much easier for :term:`FiPy` users to 
  participate in development.

.. _Git:       https://github.com/usnistgov/fipy

Fixes
-----

- `#398 <https://github.com/usnistgov/fipy/issues/398>`_:
  Home page needs out-of-NIST redirects
- `#397 <https://github.com/usnistgov/fipy/issues/397>`_:
  Switch to `sphinxcontrib-bibtex`
- `#396 <https://github.com/usnistgov/fipy/issues/396>`_:
  enable google analytics
- `#395 <https://github.com/usnistgov/fipy/issues/395>`_:
  Documentation change for Ubuntu install
- `#393 <https://github.com/usnistgov/fipy/issues/393>`_:
  `CylindricalNonUniformGrid2D` doesn't make a `FaceVariable` for
  `exteriorFaces`
- `#392 <https://github.com/usnistgov/fipy/issues/392>`_:
  `exit_nist.cgi` deprecated
- `#391 <https://github.com/usnistgov/fipy/issues/391>`_:
  Péclet inequalities have the wrong sign
- `#388 <https://github.com/usnistgov/fipy/issues/388>`_:
  Windows 64 and numpy's `dtype=int`
- `#384 <https://github.com/usnistgov/fipy/issues/384>`_:
  Add support for Matplotlib `streamplot`
- `#382 <https://github.com/usnistgov/fipy/issues/382>`_:
  Neumann boundary conditions not clearly documented
- `#381 <https://github.com/usnistgov/fipy/issues/381>`_:
  numpy 1.7.1 test failures with `physicalField.py`
- `#377 <https://github.com/usnistgov/fipy/issues/377>`_:
  `VanLeerConvectionTerm` MinMod slope limiter is broken
- `#376 <https://github.com/usnistgov/fipy/issues/376>`_:
  testing `CommitTicketUpdater`
- `#375 <https://github.com/usnistgov/fipy/issues/375>`_:
  NumPy 1.7.0 doesn't have `_formatInteger`
- `#373 <https://github.com/usnistgov/fipy/issues/373>`_:
  Bug with numpy 1.7.0
- `#372 <https://github.com/usnistgov/fipy/issues/372>`_:
  convection problem with cylindrical grid
- `#371 <https://github.com/usnistgov/fipy/issues/371>`_:
  `examples/phase/binary.py` has problems
- `#370 <https://github.com/usnistgov/fipy/issues/370>`_:
  FIPY_DISPLAY_MATRIX is broken
- `#368 <https://github.com/usnistgov/fipy/issues/368>`_:
  Viewers don't inline well in IPython notebook
- `#367 <https://github.com/usnistgov/fipy/issues/367>`_:
  Change documentation to promote use of stackoverflow
- `#366 <https://github.com/usnistgov/fipy/issues/366>`_:
  `unOps` can't be pickled
- `#365 <https://github.com/usnistgov/fipy/issues/365>`_:
  Rename communicator instances
- `#364 <https://github.com/usnistgov/fipy/issues/364>`_:
  Parallel bug in non-uniform grids and conflicting mesh class and
  factory function names
- `#360 <https://github.com/usnistgov/fipy/issues/360>`_:
  NIST CSS changed
- `#356 <https://github.com/usnistgov/fipy/issues/356>`_:
  link to mailing list is wrong
- `#353 <https://github.com/usnistgov/fipy/issues/353>`_:
  Update Ohloh to point at git repo
- `#352 <https://github.com/usnistgov/fipy/issues/352>`_:
  `getVersion()` fails on Py3k
- `#350 <https://github.com/usnistgov/fipy/issues/350>`_:
  Gmsh importer can't read mesh elements with no tags
- `#347 <https://github.com/usnistgov/fipy/issues/347>`_:
  Include mailing list activity frame on front page
- `#339 <https://github.com/usnistgov/fipy/issues/339>`_:
  Fix for test failures on `loki`
- `#337 <https://github.com/usnistgov/fipy/issues/337>`_:
  Clean up interaction between dependencies and installation process
- `#336 <https://github.com/usnistgov/fipy/issues/336>`_:
  `fipy.test()` and `fipy/test.py` clash
- `#334 <https://github.com/usnistgov/fipy/issues/334>`_:
  Make the citation links go to the DOI links
- `#333 <https://github.com/usnistgov/fipy/issues/333>`_:
  Web page links seem to be broken
- `#331 <https://github.com/usnistgov/fipy/issues/331>`_:
  Assorted errors
- `#330 <https://github.com/usnistgov/fipy/issues/330>`_:
  `faceValue` as `FaceCenters` gives inline failures
- `#329 <https://github.com/usnistgov/fipy/issues/329>`_:
  Gmsh background mesh doesn't work in parallel
- `#326 <https://github.com/usnistgov/fipy/issues/326>`_:
  `Gmsh2D` does not respect background mesh
- `#323 <https://github.com/usnistgov/fipy/issues/323>`_:
  `getFaceCenters()` should return a `FaceVariable`
- `#319 <https://github.com/usnistgov/fipy/issues/319>`_:
  Explicit convection terms should fail when the equation has no
  `TransientTerm` `(dt=None)`
- `#318 <https://github.com/usnistgov/fipy/issues/318>`_:
  FiPy will not import
- `#311 <https://github.com/usnistgov/fipy/issues/311>`_:
  LSMLIB refactor
- `#305 <https://github.com/usnistgov/fipy/issues/305>`_:
  `mpirun -np 2 python -Wd setup.py test --trilinos` hanging on
  sandbox under buildbot
- `#297 <https://github.com/usnistgov/fipy/issues/297>`_:
  Remove deprecated gist and gnuplot support
- `#291 <https://github.com/usnistgov/fipy/issues/291>`_:
  efficiency_test chokes on `liquidVapor2D.py`
- `#289 <https://github.com/usnistgov/fipy/issues/289>`_:
  `diffusionTerm._test()` requires Pysparse
- `#287 <https://github.com/usnistgov/fipy/issues/287>`_:
  move FiPy to distributed version control
- `#275 <https://github.com/usnistgov/fipy/issues/275>`_:
  `mpirun -np 2 python setup.py test --no-pysparse` hangs on `bunter`
- `#274 <https://github.com/usnistgov/fipy/issues/274>`_:
  Epetra `Norm2` failure in parallel
- `#272 <https://github.com/usnistgov/fipy/issues/272>`_:
  Error adding meshes
- `#269 <https://github.com/usnistgov/fipy/issues/269>`_:
  Rename `GridXD`
- `#255 <https://github.com/usnistgov/fipy/issues/255>`_:
  numpy 1.5.1 and masked arrays
- `#253 <https://github.com/usnistgov/fipy/issues/253>`_:
  Move the mail archive link to a more prominent place on web page.
- `#245 <https://github.com/usnistgov/fipy/issues/245>`_:
  Fix `fipy.terms._BinaryTerm` test failure in parallel
- `#228 <https://github.com/usnistgov/fipy/issues/228>`_:
  `--pysparse` configuration should never attempt MPI imports
- `#225 <https://github.com/usnistgov/fipy/issues/225>`_:
  Windows interactive plotting mostly broken
- `#209 <https://github.com/usnistgov/fipy/issues/209>`_:
  add Rhie-Chow correction term in stokes cavity example
- `#180 <https://github.com/usnistgov/fipy/issues/180>`_:
  broken arithmetic face to cell distance calculations
- `#128 <https://github.com/usnistgov/fipy/issues/128>`_:
  Trying to "solve" an integer `CellVariable` should raise an error
- `#123 <https://github.com/usnistgov/fipy/issues/123>`_:
  `numerix.dot` doesn't support tensors
- `#103 <https://github.com/usnistgov/fipy/issues/103>`_:
  `subscriber()._markStale()` `AttributeError`
- `#61 <https://github.com/usnistgov/fipy/issues/61>`_:
  Move `ImplicitDiffusionTerm().solve(var) == 0` "failure" from
  `examples.phase.simple` to `examples.diffusion.mesh1D`?

--------------------------
Version 3.0.1 - 2012-10-03
--------------------------

Fixes
-----

- `#346 <https://github.com/usnistgov/fipy/issues/346>`_:
  text in `trunk/examples/convection/source.py`
  is out of date
- `#342 <https://github.com/usnistgov/fipy/issues/342>`_:
  sign issues for equation with transient, convection and implicit
  terms
- `#338 <https://github.com/usnistgov/fipy/issues/338>`_:
  SvnToGit clean up

------------------------
Version 3.0 - 2012-08-16
------------------------

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
- FiPy is capable of running under :term:`Python 3`.
- "getter" and "setter" methods have been pervasively changed to Python 
  properties.
- The test suite now runs much faster.
- Tests can now be run on a full install using `fipy.test()`.
- The functions of the :mod:`~fipy.tools.numerix` module are no longer 
  included in the :mod:`fipy` namespace.  See
  :mod:`examples.updating.update2_0to3_0` for details.
- Equations containing a :class:`~fipy.terms.transientTerm.TransientTerm`,
  must specify the timestep by passing a ``dt=`` argument when calling
  :meth:`~fipy.terms.term.Term.solve` or :meth:`~fipy.terms.term.Term.sweep`.

.. warning::

   :term:`FiPy` 3 brought unavoidable syntax changes from :term:`FiPy` 2.
   Please see :mod:`examples.updating.update2_0to3_0` for guidance on the
   changes that you will need to make to your :term:`FiPy` 2.x scripts.

Fixes
-----

- `#332 <https://github.com/usnistgov/fipy/issues/332>`_:
  Inline failure on Ubuntu x86_64
- `#324 <https://github.com/usnistgov/fipy/issues/324>`_:
  constraining values with `ImplicitSourceTerm` not documented?
- `#317 <https://github.com/usnistgov/fipy/issues/317>`_:
  `gmshImport` tests fail on Windows due to shared file
- `#316 <https://github.com/usnistgov/fipy/issues/316>`_:
  changes to `gmshImport.py` caused `--inline` problems
- `#313 <https://github.com/usnistgov/fipy/issues/313>`_:
  Gmsh I/O
- `#307 <https://github.com/usnistgov/fipy/issues/307>`_:
  Failures on sandbox under buildbot
- `#306 <https://github.com/usnistgov/fipy/issues/306>`_:
  Add in parallel buildbot testing on more than 2 processors
- `#302 <https://github.com/usnistgov/fipy/issues/302>`_:
  `CellVariable.min()` broken in parallel
- `#301 <https://github.com/usnistgov/fipy/issues/301>`_:
  `Epetra.PyComm()` broken on Debian
- `#300 <https://github.com/usnistgov/fipy/issues/300>`_:
  `examples/cahnHilliard/mesh2D.py` broken with -- trilinos
- `#299 <https://github.com/usnistgov/fipy/issues/299>`_:
  Viewers not working when plotting meshes with zero cells in parallel
- `#298 <https://github.com/usnistgov/fipy/issues/298>`_:
  Memory consumption growth with repeated meshing, especially with
  Gmsh
- `#294 <https://github.com/usnistgov/fipy/issues/294>`_:
  `--pysparse --inline` failures
- `#293 <https://github.com/usnistgov/fipy/issues/293>`_:
  `python examples/cahnHilliard/sphere.py --inline` segfaults on OS X
- `#292 <https://github.com/usnistgov/fipy/issues/292>`_:
  two `--scipy` failures
- `#290 <https://github.com/usnistgov/fipy/issues/290>`_:
  Improve test reporting to avoid inconsequential buildbot failures
- `#288 <https://github.com/usnistgov/fipy/issues/288>`_:
  gmsh importer and gmsh tests don't clean up after themselves
- `#286 <https://github.com/usnistgov/fipy/issues/286>`_:
  get running in Py3k
- `#285 <https://github.com/usnistgov/fipy/issues/285>`_:
  remove deprecated `viewers.make()`
- `#284 <https://github.com/usnistgov/fipy/issues/284>`_:
  remove deprecated `Variable.transpose()`
- `#281 <https://github.com/usnistgov/fipy/issues/281>`_:
  remove deprecated `NthOrderDiffusionTerm`
- `#280 <https://github.com/usnistgov/fipy/issues/280>`_:
  remove deprecated `diffusionTerm=` argument to `ConvectionTerm`
- `#277 <https://github.com/usnistgov/fipy/issues/277>`_:
  remove deprecated `steps=` from Solver
- `#273 <https://github.com/usnistgov/fipy/issues/273>`_:
  Make `DiffusionTermNoCorrection` the default
- `#270 <https://github.com/usnistgov/fipy/issues/270>`_:
  tests take *too* long!!!
- `#267 <https://github.com/usnistgov/fipy/issues/267>`_:
  Reduce the run times for chemotaxis tests
- `#264 <https://github.com/usnistgov/fipy/issues/264>`_:
  HANG in parallel test of `examples/chemotaxis/input2D.py` on some
  configurations
- `#261 <https://github.com/usnistgov/fipy/issues/261>`_:
  `GmshImport` should read element colors
- `#260 <https://github.com/usnistgov/fipy/issues/260>`_:
  `GmshImport` should support all element types
- `#259 <https://github.com/usnistgov/fipy/issues/259>`_:
  Introduce `mesh.x` as shorthand for `mesh.cellCenters[0]` etc
- `#258 <https://github.com/usnistgov/fipy/issues/258>`_:
  `GmshExport` is not tested and does not work
- `#252 <https://github.com/usnistgov/fipy/issues/252>`_:
  Include Benny's improved interpolation patch
- `#250 <https://github.com/usnistgov/fipy/issues/250>`_:
  TeX is wrong in `examples.phase.quaternary`
- `#247 <https://github.com/usnistgov/fipy/issues/247>`_:
  `diffusionTerm(var=var1).solver(var=var0)` should fail sensibly
- `#243 <https://github.com/usnistgov/fipy/issues/243>`_:
  close out reconstrain branch
- `#242 <https://github.com/usnistgov/fipy/issues/242>`_:
  update documentation
- `#240 <https://github.com/usnistgov/fipy/issues/240>`_:
  Profile and merge reconstrain branch
- `#237 <https://github.com/usnistgov/fipy/issues/237>`_:
  `--Trilinos --no-pysparse` uses Pysparse?!?
- `#236 <https://github.com/usnistgov/fipy/issues/236>`_:
  anisotropic diffusion and constraints don't mix
- `#235 <https://github.com/usnistgov/fipy/issues/235>`_:
  changed constraints don't propagate
- `#231 <https://github.com/usnistgov/fipy/issues/231>`_:
  `factoryMeshes.py` not up to date with respect to keyword arguments
- `#223 <https://github.com/usnistgov/fipy/issues/223>`_:
  mesh in FiPy name space
- `#218 <https://github.com/usnistgov/fipy/issues/218>`_:
  Absence of `enthought.tvtk` causes test failures
- `#216 <https://github.com/usnistgov/fipy/issues/216>`_:
  Fresh FiPy gives "`ImportError: No viewers found`"
- `#213 <https://github.com/usnistgov/fipy/issues/213>`_:
  PyPI is failing
- `#206 <https://github.com/usnistgov/fipy/issues/206>`_:
  `gnuplot1d` gives error on plot of `FaceVariable`
- `#205 <https://github.com/usnistgov/fipy/issues/205>`_:
  wrong cell to cell normal in periodic meshes
- `#203 <https://github.com/usnistgov/fipy/issues/203>`_:
  Give helpful error on - or / of meshes
- `#202 <https://github.com/usnistgov/fipy/issues/202>`_:
  mesh manipulation of periodic meshes leads to errors
- `#201 <https://github.com/usnistgov/fipy/issues/201>`_:
  Use physical velocity in the manual/FAQ
- `#200 <https://github.com/usnistgov/fipy/issues/200>`_:
  FAQ gives bad guidance for anisotropic diffusion
- `#195 <https://github.com/usnistgov/fipy/issues/195>`_:
  term multiplication changes result
- `#163 <https://github.com/usnistgov/fipy/issues/163>`_:
  Default time steps should be infinite
- `#162 <https://github.com/usnistgov/fipy/issues/162>`_:
  remove ones and zeros from `numerix.py`
- `#130 <https://github.com/usnistgov/fipy/issues/130>`_:
  tests should be run with `fipy.tests()`
- `#86 <https://github.com/usnistgov/fipy/issues/86>`_:
  Grids should take `Lx`, `Ly`, `Lz` arguments
- `#77 <https://github.com/usnistgov/fipy/issues/77>`_:
  `CellVariable.hasOld()` should set `self.old`
- `#44 <https://github.com/usnistgov/fipy/issues/44>`_:
  Navier-Stokes

--------------------------
Version 2.1.3 - 2012-01-17
--------------------------

Fixes
-----

- `#282 <https://github.com/usnistgov/fipy/issues/282>`_:
  remove deprecated getters and setters
- `#279 <https://github.com/usnistgov/fipy/issues/279>`_:
  remove deprecated `fipy.meshes.numMesh` submodule
- `#278 <https://github.com/usnistgov/fipy/issues/278>`_:
  remove deprecated forms of Gmsh meshes
- `#268 <https://github.com/usnistgov/fipy/issues/268>`_:
  Set up `Zizou` as a working slave
- `#262 <https://github.com/usnistgov/fipy/issues/262>`_:
  issue with solvers
- `#256 <https://github.com/usnistgov/fipy/issues/256>`_:
  `Grid1D(dx=(1,2,3))` failure
- `#251 <https://github.com/usnistgov/fipy/issues/251>`_:
  parallel is broken
- `#241 <https://github.com/usnistgov/fipy/issues/241>`_:
  Set Sandbox up as a working slave
- `#238 <https://github.com/usnistgov/fipy/issues/238>`_:
  `_BinaryTerm.var` is not predictable
- `#233 <https://github.com/usnistgov/fipy/issues/233>`_:
  coupled convection-diffusion always treated as Upwind
- `#224 <https://github.com/usnistgov/fipy/issues/224>`_:
  "matrices are not aligned" errors in example test suite
- `#222 <https://github.com/usnistgov/fipy/issues/222>`_:
  Non-uniform `Grid3D` fails to __add__
- `#221 <https://github.com/usnistgov/fipy/issues/221>`_:
  Problem with fipy and gmsh
- `#219 <https://github.com/usnistgov/fipy/issues/219>`_:
  matforge css is hammer-headed
- `#208 <https://github.com/usnistgov/fipy/issues/208>`_:
  numpy 2.0: `arrays have a dot method`
- `#207 <https://github.com/usnistgov/fipy/issues/207>`_:
  numpy 2.0: `masked arrays cast right of product to ndarray`
- `#196 <https://github.com/usnistgov/fipy/issues/196>`_:
  Pysparse won't import in Python 2.6.5 on Windows
- `#152 <https://github.com/usnistgov/fipy/issues/152>`_:
  (Re)Implement SciPy solvers
- `#138 <https://github.com/usnistgov/fipy/issues/138>`_:
  FAQ on boundary conditions
- `#100 <https://github.com/usnistgov/fipy/issues/100>`_:
  testing from the Windows dist using the ipython command line
- `#80 <https://github.com/usnistgov/fipy/issues/80>`_:
  Windows - testing - idle `-ipython`
- `#46 <https://github.com/usnistgov/fipy/issues/46>`_:
  Variable needs to consider boundary conditions
- `#45 <https://github.com/usnistgov/fipy/issues/45>`_:
  Slicing a vector Variable should produce a scalar Variable

--------------------------
Version 2.1.2 - 2011-04-20
--------------------------

The significant changes since version 2.1.1 are:

- :term:`Trilinos` efficiency improvements
- Diagnostics of the parallel environment

Fixes
-----

- `#232 <https://github.com/usnistgov/fipy/issues/232>`_:
  Mayavi broken on windows because it has no `SIGHUP`.
- `#230 <https://github.com/usnistgov/fipy/issues/230>`_:
  `factoryMeshes.py` not up to date with respect to keyword arguments
- `#226 <https://github.com/usnistgov/fipy/issues/226>`_:
  `MatplotlibViewer` fails if backend doesn't support `flush_events()`
- `#225 <https://github.com/usnistgov/fipy/issues/225>`_:
  Windows interactive plotting mostly broken
- `#217 <https://github.com/usnistgov/fipy/issues/217>`_:
  Gmsh `CellVariables` can't be unpickled
- `#191 <https://github.com/usnistgov/fipy/issues/191>`_:
  `sphereDaemon.py` missing in FiPy 2.1 and from trunk
- `#187 <https://github.com/usnistgov/fipy/issues/187>`_:
  Concatenated `Mesh` garbled by `dump.write`/`read`

--------------------------
Version 2.1.1 - 2010-10-05
--------------------------

The significant changes since version 2.1 are:

- :class:`~fipy.viewers.matplotlibViewer.MatplotlibViewer` can display 
  into an existing set of Matplotlib axes.

- :term:`Pysparse` and :term:`Trilinos` are now completely independent.

Fixes
-----

- `#199 <https://github.com/usnistgov/fipy/issues/199>`_:
  dummy viewer results in
  "`NotImplementedError: can't instantiate abstract base class`"
- `#198 <https://github.com/usnistgov/fipy/issues/198>`_:
  bug problem with `CylindricalGrid1D`
- `#197 <https://github.com/usnistgov/fipy/issues/197>`_:
  How to tell if parallel is configured properly?
- `#194 <https://github.com/usnistgov/fipy/issues/194>`_:
  `FIPY_DISPLAY_MATRIX` on empty matrix with large b-vector throws
  `ValueError`
- `#193 <https://github.com/usnistgov/fipy/issues/193>`_:
  `FIPY_DISPLAY_MATRIX` raises `ImportError` in FiPy 2.1 and trunk
- `#192 <https://github.com/usnistgov/fipy/issues/192>`_:
  `FIPY_DISPLAY_MATRIX=terms` raises `TypeError` in FiPy 2.1 and trunk

------------------------
Version 2.1 - 2010-04-01
------------------------

The relatively small change in version number belies significant advances
in :term:`FiPy` capabilities.  This release did not receive a "full"
version increment because it is completely (er...  [#almost]_) compatible
with older scripts.

The significant changes since version 2.0.2 are:

- :term:`FiPy` can use :term:`Trilinos` for :ref:`PARALLEL`.

- We have switched from :term:`MayaVi` 1 to :term:`Mayavi` 2. This 
  :class:`~fipy.viewers.viewer.Viewer` is an independent process that 
  allows interaction with the display while a simulation is running.

- Documentation has been switched to :term:`Sphinx`, allowing the entire manual
  to be available on the web and for our documentation to link to the
  documentation for packages such as :mod:`numpy`, :mod:`scipy`,
  :mod:`matplotlib`, and for :term:`Python` itself.

Fixes
-----

- `#190 <https://github.com/usnistgov/fipy/issues/190>`_:
  "matplotlib: list index out of range" when no title given, but only
  sometimes
- `#182 <https://github.com/usnistgov/fipy/issues/182>`_:
  `~binOp` doesn't work on branches/version-2_0
- `#180 <https://github.com/usnistgov/fipy/issues/180>`_:
  broken arithmetic face to cell distance calculations
- `#179 <https://github.com/usnistgov/fipy/issues/179>`_:
  `easy_install` instructions for Mac OS X are broken
- `#177 <https://github.com/usnistgov/fipy/issues/177>`_:
  broken `setuptools` url with python 2.6
- `#169 <https://github.com/usnistgov/fipy/issues/169>`_:
  The FiPy webpage seems to be broken on Internet Explorer
- `#156 <https://github.com/usnistgov/fipy/issues/156>`_:
  update the mayavi viewer to use  mayavi 2
- `#153 <https://github.com/usnistgov/fipy/issues/153>`_:
  Switch documentation to use `:math:` directive

.. [#almost] Only two examples from :term:`FiPy` 2.0 fail when run with
   :term:`FiPy` 2.1:

    * :file:`examples/phase/symmetry.py` fails because
      :class:`~fipy.meshes.mesh.Mesh` no longer provides a
      ``getCells`` method. The mechanism
      for enforcing symmetry in the updated example is both clearer and
      faster.

    * :mod:`examples.levelSet.distanceFunction.circle` fails because of a
      change in the comparison of masked values.

   Both of these are subtle issues unlikely to affect very many
   :term:`FiPy` users.

--------------------------
Version 2.0.3 - 2010-03-17
--------------------------

Fixes
-----

- `#188 <https://github.com/usnistgov/fipy/issues/188>`_:
  `SMTPSenderRefused: (553, "5.1.8 <trac@matdl-osi.org>... Domain of sender address trac@matdl-osi.org does not exist", u'"FiPy" <trac@matdl-osi.org>')`
- `#184 <https://github.com/usnistgov/fipy/issues/184>`_:
  `gmshExport.exportAsMesh()` doesn't work
- `#183 <https://github.com/usnistgov/fipy/issues/183>`_:
  FiPy 2.0.2 `LinearJORSolver.__init__`  calls `Solver` rather than
  `PysparseSolver`
- `#181 <https://github.com/usnistgov/fipy/issues/181>`_:
  Navier-Stokes again
- `#151 <https://github.com/usnistgov/fipy/issues/151>`_:
  update mayavi viewer to use mayavi2
- `#13 <https://github.com/usnistgov/fipy/issues/13>`_:
  Mesh refactor

--------------------------
Version 2.0.2 - 2009-06-11
--------------------------

Fixes
-----

- `#176 <https://github.com/usnistgov/fipy/issues/176>`_:
  Win32 distribution test error
- `#175 <https://github.com/usnistgov/fipy/issues/175>`_:
  `Grid3D` `getFaceCenters` incorrect when mesh is offset
- `#170 <https://github.com/usnistgov/fipy/issues/170>`_:
  `Variable` doesn't implement `__invert__`

--------------------------
Version 2.0.1 - 2009-04-23
--------------------------

Fixes
-----

- `#154 <https://github.com/usnistgov/fipy/issues/154>`_:
  Update manuals

------------------------
Version 2.0 - 2009-02-09
------------------------

.. warning::

   :term:`FiPy` 2 brings unavoidable syntax changes. Please see
   :mod:`examples.updating.update1_0to2_0` for guidance on the changes that
   you will need to make to your :term:`FiPy` 1.x scripts.

The significant changes since version 1.2 are:

- :class:`~fipy.variables.cellVariable.CellVariable` and
  :class:`~fipy.variables.faceVariable.FaceVariable` objects can hold
  values of any rank.

- Much simpler syntax for specifying
  ``Cell``\s for initial conditions and
  ``Face``\s for boundary conditions.

- Automated determination of the Péclet number and partitioning of 
  :class:`~fipy.terms.implicitSourceTerm.ImplicitSourceTerm` coefficients
  between the matrix diagonal and the right-hand-side-vector.

- Simplified :class:`~fipy.viewers.viewer.Viewer` syntax.

- Support for the `Trilinos solvers`_.

- Support for anisotropic diffusion coefficients.

.. _Trilinos solvers: http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://trilinos.sandia.gov

- `#167 <https://github.com/usnistgov/fipy/issues/167>`_:
  example showing how to go from 1.2 to 2.0
- `#166 <https://github.com/usnistgov/fipy/issues/166>`_:
  Still references to `VectorCell` and `VectorFace` `Variable` in manual
- `#165 <https://github.com/usnistgov/fipy/issues/165>`_:
  Edit the what's new section of the manual
- `#149 <https://github.com/usnistgov/fipy/issues/149>`_:
  Test viewers
- `#143 <https://github.com/usnistgov/fipy/issues/143>`_:
  Document syntax changes
- `#141 <https://github.com/usnistgov/fipy/issues/141>`_:
  enthought toolset?
- `#140 <https://github.com/usnistgov/fipy/issues/140>`_:
  easy_install fipy
- `#136 <https://github.com/usnistgov/fipy/issues/136>`_:
  Document anisotropic diffusion
- `#135 <https://github.com/usnistgov/fipy/issues/135>`_:
  Trilinos documentation
- `#127 <https://github.com/usnistgov/fipy/issues/127>`_:
  Examples can be very fragile with respect to floating point

-------------------------
Version 1.2.3 - 2009-01-0
-------------------------

Fixes
-----

- `#54 <https://github.com/usnistgov/fipy/issues/54>`_:
  `python setup.py test` fails

--------------------------
Version 1.2.2 - 2008-12-30
--------------------------

Fixes
-----

- `#161 <https://github.com/usnistgov/fipy/issues/161>`_:
  get pysparse working with python 2.4
- `#160 <https://github.com/usnistgov/fipy/issues/160>`_:
  Grid class
- `#157 <https://github.com/usnistgov/fipy/issues/157>`_:
  temp files on widows
- `#155 <https://github.com/usnistgov/fipy/issues/155>`_:
  fix some of the deprecation warnings appearing in the tests
- `#150 <https://github.com/usnistgov/fipy/issues/150>`_:
  PythonXY installation?
- `#148 <https://github.com/usnistgov/fipy/issues/148>`_:
  SciPy 0.7.0 solver failures on Macs
- `#147 <https://github.com/usnistgov/fipy/issues/147>`_:
  Disable CGS solver in pysparse
- `#145 <https://github.com/usnistgov/fipy/issues/145>`_:
  `Viewer` factory fails for rank-1 `CellVariable`
- `#144 <https://github.com/usnistgov/fipy/issues/144>`_:
  intermittent failure on 
  `examples/diffusion/explicit/mixedelement.py --inline`
- `#142 <https://github.com/usnistgov/fipy/issues/142>`_:
  merge Viewers branch
- `#139 <https://github.com/usnistgov/fipy/issues/139>`_:
  Get a Windows Bitten build slave
- `#137 <https://github.com/usnistgov/fipy/issues/137>`_:
  Backport examples from manuscript
- `#131 <https://github.com/usnistgov/fipy/issues/131>`_:
  `MatplotlibViewer` doesn't properly report the supported file
  extensions
- `#126 <https://github.com/usnistgov/fipy/issues/126>`_:
  Variable, float, integer
- `#125 <https://github.com/usnistgov/fipy/issues/125>`_:
  Pickled test data embeds obsolete packages
- `#124 <https://github.com/usnistgov/fipy/issues/124>`_:
  Can't pickle a `binOp`
- `#121 <https://github.com/usnistgov/fipy/issues/121>`_:
  `simpleTrenchSystem.py`
- `#120 <https://github.com/usnistgov/fipy/issues/120>`_:
  mayavi display problems
- `#118 <https://github.com/usnistgov/fipy/issues/118>`_:
  Automatically handle casting of `Variable` from `int` to `float`
  when necessary.
- `#117 <https://github.com/usnistgov/fipy/issues/117>`_:
  `getFacesBottom`, `getFacesTop` etc. lack clear description in the
  reference
- `#115 <https://github.com/usnistgov/fipy/issues/115>`_:
  viewing 3D Cahn-Hilliard is broken
- `#113 <https://github.com/usnistgov/fipy/issues/113>`_:
  OS X (MacBook Pro; Intel) FiPy installation problems
- `#112 <https://github.com/usnistgov/fipy/issues/112>`_:
  `stokesCavity.py` doesn't display properly with matplotlib
- `#111 <https://github.com/usnistgov/fipy/issues/111>`_:
  Can't display `Grid2D` variables with matplotlib on Linux
- `#110 <https://github.com/usnistgov/fipy/issues/110>`_:
  "Numeric array value must be dimensionless"  in ElPhF examples
- `#109 <https://github.com/usnistgov/fipy/issues/109>`_:
  doctest of `fipy.variables.variable.Variable.__array__`
- `#108 <https://github.com/usnistgov/fipy/issues/108>`_:
  `numerix.array * FaceVariable` is broken
- `#107 <https://github.com/usnistgov/fipy/issues/107>`_:
  Can't move matplotlib windows on Mac
- `#106 <https://github.com/usnistgov/fipy/issues/106>`_:
  Concatenation of `Grid1D` objects doesn't always work
- `#105 <https://github.com/usnistgov/fipy/issues/105>`_:
  useless broken __array__ tests should be removed
- `#102 <https://github.com/usnistgov/fipy/issues/102>`_:
  viewer limits should just be set as arguments, rather than as a dict
- `#99 <https://github.com/usnistgov/fipy/issues/99>`_:
  `Matplotlib2DGridViewer` cannot update multiple views
- `#97 <https://github.com/usnistgov/fipy/issues/97>`_:
  Windows does not seem to handle NaN correctly.
- `#96 <https://github.com/usnistgov/fipy/issues/96>`_:
  broken tests with version 2.0 of gmsh
- `#95 <https://github.com/usnistgov/fipy/issues/95>`_:
  attached code breaks with `--inline`
- `#92 <https://github.com/usnistgov/fipy/issues/92>`_:
  Pygist is dead (it's official)
- `#84 <https://github.com/usnistgov/fipy/issues/84>`_:
  Test failures on Intel Mac
- `#83 <https://github.com/usnistgov/fipy/issues/83>`_:
  `ZeroDivisionError` for `CellTerm` when calling `getOld()` on its
  coefficient
- `#79 <https://github.com/usnistgov/fipy/issues/79>`_:
  `viewers.make()` to `viewers.Viewer()`
- `#67 <https://github.com/usnistgov/fipy/issues/67>`_:
  Mesh viewing and unstructured data.
- `#43 <https://github.com/usnistgov/fipy/issues/43>`_:
  `TSVViewer` doesn't always get the right shape for the var
- `#34 <https://github.com/usnistgov/fipy/issues/34>`_:
  float(&infinity&) issue on windows

--------------------------
Version 1.2.1 - 2008-02-08
--------------------------

Fixes
-----

- `#122 <https://github.com/usnistgov/fipy/issues/122>`_:
  check argument types for meshes
- `#119 <https://github.com/usnistgov/fipy/issues/119>`_:
  max is broken for Variables
- `#116 <https://github.com/usnistgov/fipy/issues/116>`_:
  Linux: failed test, `TypeError: No array interface...` in solve()
- `#104 <https://github.com/usnistgov/fipy/issues/104>`_:
  Syntax error in `MatplotlibVectorViewer._plot()`
- `#101 <https://github.com/usnistgov/fipy/issues/101>`_:
  matplotlib 1D viewer autoscales when a limit is set to 0
- `#93 <https://github.com/usnistgov/fipy/issues/93>`_:
  Broken examples
- `#91 <https://github.com/usnistgov/fipy/issues/91>`_:
  update the examples to use `from fipy import *`
- `#76 <https://github.com/usnistgov/fipy/issues/76>`_:
  `solve()` and `sweep()` accept `dt=CellVariable`
- `#75 <https://github.com/usnistgov/fipy/issues/75>`_:
  installation of fipy should auto include README as a docstring
- `#74 <https://github.com/usnistgov/fipy/issues/74>`_:
  Some combinations of `DiffusionTerm` and `ConvectionTerm` do not work
- `#51 <https://github.com/usnistgov/fipy/issues/51>`_:
  __pos__ doesn't work for terms
- `#50 <https://github.com/usnistgov/fipy/issues/50>`_:
  Broken examples
- `#39 <https://github.com/usnistgov/fipy/issues/39>`_:
  matplotlib broken on mac with version 0.72.1
- `#19 <https://github.com/usnistgov/fipy/issues/19>`_:
  Péclet number
- `#15 <https://github.com/usnistgov/fipy/issues/15>`_:
  Boundary conditions and Terms

------------------------
Version 1.2 - 2007-02-12
------------------------

The significant changes since version 1.1 are:

- `--inline` automatically generates C code from `Variable` expressions.

- :term:`FiPy` has been updated to use the :term:`Python` :term:`NumPy` module.
  :term:`FiPy` no longer works with the older :term:`Numeric` module.

Fixes
-----

- `#98 <https://github.com/usnistgov/fipy/issues/98>`_:
  Windows patch for some broken test cases
- `#94 <https://github.com/usnistgov/fipy/issues/94>`_:
  `--inline` error for attached code
- `#90 <https://github.com/usnistgov/fipy/issues/90>`_:
  bug in matplotlib 0.87.7:
  `TypeError: only length-1 arrays can be converted to Python scalars`.
- `#72 <https://github.com/usnistgov/fipy/issues/72>`_:
  needless rebuilding of variables
- `#66 <https://github.com/usnistgov/fipy/issues/66>`_:
  PDF rendering issues for the guide on various platforms
- `#62 <https://github.com/usnistgov/fipy/issues/62>`_:
  fipy guide pdf bug: "`an unrecognized token 13c was found`"
- `#55 <https://github.com/usnistgov/fipy/issues/55>`_:
  Error for internal BCs
- `#52 <https://github.com/usnistgov/fipy/issues/52>`_:
  `FaceVariable * FaceVectorVariable` memory
- `#48 <https://github.com/usnistgov/fipy/issues/48>`_:
  Documentation is not inherited from &hidden& classes
- `#42 <https://github.com/usnistgov/fipy/issues/42>`_:
  `fipy.models.phase.phase.addOverFacesVariable` is gross
- `#41 <https://github.com/usnistgov/fipy/issues/41>`_:
  :file:`EFFICIENCY.txt` example fails to make viewer
- `#30 <https://github.com/usnistgov/fipy/issues/30>`_:
  periodic boundary condition support
- `#25 <https://github.com/usnistgov/fipy/issues/25>`_:
  make phase field examples more explicit
- `#23 <https://github.com/usnistgov/fipy/issues/23>`_:
  sweep control, iterator object, error norms
- `#21 <https://github.com/usnistgov/fipy/issues/21>`_:
  Update FiPy to use numpy
- `#16 <https://github.com/usnistgov/fipy/issues/16>`_:
  Dimensions
- `#12 <https://github.com/usnistgov/fipy/issues/12>`_:
  Refactor viewers
- `#1 <https://github.com/usnistgov/fipy/issues/1>`_:
  Gnuplot doesn't display on windows

------------------------
Version 1.1 - 2006-06-06
------------------------

The significant changes since version 1.0 are:

- Memory efficiency has been improved in a number of ways, but most
  significantly by:

  * not caching all intermediate ``Variable`` values.
  * introducing ``UniformGrid`` classes that calculate geometric
    arrays on the fly.

  Details of these improvements are presented in :ref:`chap:Efficiency`.

- Installation on Windows has been made considerably easier by
  constructing executable installers for :term:`FiPy` and its
  dependencies.

- The arithmetic for ``Variable`` subclasses now works, and returns
  sensible answers. For example, ``VectorCellVariable * CellVariable``
  returns a ``VectorCellVariable``.

- ``PeriodicGrid`` meshes have been implemented. Currently, however,
  there and no examples of their use in the manual.

- Many of the examples have been completely rewritten

  * A basic 1D diffusion problem now serves as a general tutorial for 
    setting up any problem in :term:`FiPy`. 
  * Several more phase field examples have been added that should make it 
    clearer how to get from the simple 1D case to the more elaborate 
    multicomponent, multidimensional, and anisotropic models.
  * The "Superfill" examples have been substantially improved with better
    functionality and documentation.
  * An example of fluid flow with the classic Stokes moving lid has been 
    added.

- A clear distinction has been made between solving an equation via `solve()`
  and iterating an non-linear equation to solution via `sweep()`. An extensive 
  explanation of the concepts involved has been added to the :ref:`FAQ`.

- Added a `MultiViewer` class that automatically groups several viewers 
  together if the variables couldn't be displayed by a single viewer.

- The abbreviated syntax ``from fipy import Class`` or ``from fipy import *``
  promised in version 1.0 actually works now. The examples all still use the
  fully qualified names.

- The repository has been converted from a CVS to a Subversion_
  repository. Details on how to check out the new repository are given
  in :ref:`INSTALLATION`.

- The :term:`FiPy` repository has also been moved from Sourceforge_ to the
  `Materials Digital Library Pathway`_.

..  _Subversion: https://subversion.apache.org/
..  _Sourceforge: https://sourceforge.net/
..  _Materials Digital Library Pathway: https://www.kent.edu/cmi/materials-digital-library-pathway-matdl

------------------------
Version 1.0 - 2005-09-16
------------------------

Numerous changes have been made since :term:`FiPy` 0.1 was released, but the most
significant ones are:

- ``Equation`` objects no longer exist. PDEs are constructed from ``Term`` 
  objects. ``Term`` objects can be added, subtracted, and equated to build up 
  an equation.

- A true 1D grid class has been added: ``fipy.meshes.grid1D.Grid1D``.

- A generic "factory" method ``fipy.viewers.make()`` has been added that will 
  do a reasonable job of automatically creating a ``Viewer`` for the supplied 
  ``Variable`` objects. The ``FIPY_VIEWER`` environment variable allows you to 
  specify your preferred viewer.

- A simple ``TSVViewer`` has been added to allow display or export to a file of 
  your solution data.

- It is no longer necessary to ``transpose()`` scalar fields in order to 
  multiply them with vector fields.

- Better default choice of solver when convection is present.

- Better examples.

- A number of `NoiseVariable` objects have been added.

- A new viewer based on :term:`Matplotlib` has been added.

- The `PyX` viewer has been removed.

- Considerably simplified the public interface to FiPy.

- Support for Python 2.4.

- Improved layout of the manuals.

- ``getLaplacian()`` method has been removed from ``CellVariable`` objects.
  You can obtain the same effect with ``getFaceGrad().getDivergence()``, 
  which provides better control.

- An ``import`` shorthand has been added that allows for::

     from fipy import Class

  instead of::

     from fipy.some.deeply.nested.module.class import Class

  This system is still experimental. Please tell us if you find situations
  that don't work.

The syntax of :term:`FiPy` 1.0 scripts is incompatible with earlier
releases.  A tutorial for updating your existing scripts can be found in
:file:`examples/updating/update0_1to1_0.py`.

Fixes
-----

- `#49 <https://github.com/usnistgov/fipy/issues/49>`_:
  Documentation for many `ConvectionTerms` is wrong
- `#47 <https://github.com/usnistgov/fipy/issues/47>`_:
  Terms should throw an error on bad `coeff` type
- `#40 <https://github.com/usnistgov/fipy/issues/40>`_:
  broken levelset test case
- `#38 <https://github.com/usnistgov/fipy/issues/38>`_:
  multiple BCs on one face broken?
- `#37 <https://github.com/usnistgov/fipy/issues/37>`_:
  Better support for periodic boundary conditions
- `#36 <https://github.com/usnistgov/fipy/issues/36>`_:
  Gnuplot doesn't display the :file:`~examples/levelSet/electroChem`
  problem on windows.
- `#35 <https://github.com/usnistgov/fipy/issues/35>`_:
  gmsh write problem on windows
- `#33 <https://github.com/usnistgov/fipy/issues/33>`_:
  `DiffusionTerm(coeff = CellVariable)` functionality
- `#32 <https://github.com/usnistgov/fipy/issues/32>`_:
  conflict_handler = `ignore` not valid in Python 2.4
- `#31 <https://github.com/usnistgov/fipy/issues/31>`_:
  Support simple import notation
- `#29 <https://github.com/usnistgov/fipy/issues/29>`_:
  periodic boundary conditions are broken
- `#28 <https://github.com/usnistgov/fipy/issues/28>`_:
  invoke the == for terms
- `#26 <https://github.com/usnistgov/fipy/issues/26>`_:
  doctest extraction with python 2.4
- `#24 <https://github.com/usnistgov/fipy/issues/24>`_:
  Pysparse windows binaries
- `#22 <https://github.com/usnistgov/fipy/issues/22>`_:
  automated efficiency_test problems
- `#20 <https://github.com/usnistgov/fipy/issues/20>`_:
  Test with Python version 2.4
- `#18 <https://github.com/usnistgov/fipy/issues/18>`_:
  Memory leak for the leveling problem
- `#17 <https://github.com/usnistgov/fipy/issues/17>`_:
  `distanceVariable` is broken
- `#14 <https://github.com/usnistgov/fipy/issues/14>`_:
  Testing mailing list interface
- `#11 <https://github.com/usnistgov/fipy/issues/11>`_:
  Reconcile versions of pysparse
- `#10 <https://github.com/usnistgov/fipy/issues/10>`_:
  check phase field crystal growth
- `#9 <https://github.com/usnistgov/fipy/issues/9>`_:
  implement levelling surfactant equation
- `#8 <https://github.com/usnistgov/fipy/issues/8>`_:
  merge `depositionRateVar` and `extensionVelocity`
- `#7 <https://github.com/usnistgov/fipy/issues/7>`_:
  Automate FiPy efficiency test
- `#6 <https://github.com/usnistgov/fipy/issues/6>`_:
  FiPy breaks on windows with Numeric 23.6
- `#5 <https://github.com/usnistgov/fipy/issues/5>`_:
  axisymmetric 2D mesh
- `#4 <https://github.com/usnistgov/fipy/issues/4>`_:
  Windows installation wizard
- `#3 <https://github.com/usnistgov/fipy/issues/3>`_:
  Windows installation instructions
- `#2 <https://github.com/usnistgov/fipy/issues/2>`_:
  Some tests fail on windows XP

-------------
Version 0.1.1
-------------

------------------------
Version 0.1 - 2004-11-05
------------------------

Original release
