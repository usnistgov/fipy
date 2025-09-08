.. |.continuousintegration| replace:: Continuous Integration
.. _.continuousintegration: https://pages.nist.gov/fipy/en/latest/ADMINISTRATA.html#continuousintegration
.. |.create_conda_environment| replace:: Create a conda_ environment
.. _.create_conda_environment: https://pages.nist.gov/fipy/en/latest/INSTALLATION.html#create-conda-environment
.. |.documentation-colon-git| replace:: Git usage
.. _.documentation-colon-git: https://pages.nist.gov/fipy/en/latest/INSTALLATION.html#documentation-git
.. |.environment| replace:: Development Environment
.. _.environment: https://pages.nist.gov/fipy/en/latest/INSTALLATION.html#environment
.. |.FiPy| replace:: FiPy
.. _.FiPy: https://pages.nist.gov/fipy/en/latest/glossary.html#term-FiPy
.. |.flagsandenvironmentvariables| replace:: Command-line Flags and Environment Variables
.. _.flagsandenvironmentvariables: https://pages.nist.gov/fipy/en/latest/USAGE.html#flagsandenvironmentvariables
.. |.Gmsh| replace:: Gmsh
.. _.Gmsh: https://pages.nist.gov/fipy/en/latest/glossary.html#term-Gmsh
.. |.installation| replace:: Installation
.. _.installation: https://pages.nist.gov/fipy/en/latest/INSTALLATION.html#installation
.. |.NumPy| replace:: NumPy
.. _.NumPy: https://pages.nist.gov/fipy/en/latest/glossary.html#term-NumPy
.. |.optionalpackages| replace:: Optional Packages
.. _.optionalpackages: https://pages.nist.gov/fipy/en/latest/INSTALLATION.html#optionalpackages
.. |.part-colon-examples+examples| replace:: examples
.. _.part-colon-examples+examples: https://pages.nist.gov/fipy/en/latest/EXAMPLES.html#part-examples
.. |.pip| replace:: pip
.. _.pip: https://pages.nist.gov/fipy/en/latest/glossary.html#term-pip
.. |.Python| replace:: Python
.. _.Python: https://pages.nist.gov/fipy/en/latest/glossary.html#term-Python
.. |.SciPy| replace:: SciPy
.. _.SciPy: https://pages.nist.gov/fipy/en/latest/glossary.html#term-SciPy
.. |.solvers| replace:: Solvers
.. _.solvers: https://pages.nist.gov/fipy/en/latest/SOLVERS.html#solvers
.. |.usage| replace:: Using FiPy
.. _.usage: https://pages.nist.gov/fipy/en/latest/USAGE.html#usage
.. |.viewers| replace:: Viewers
.. _.viewers: https://pages.nist.gov/fipy/en/latest/VIEWERS.html#viewers


.. _INSTALLATION:

============
Installation
============

The |.FiPy|_ finite volume PDE solver relies on several
third-party packages.  It is *best to obtain and install those first*
before attempting to install |.FiPy|_. This document explains how
to install |.FiPy|_, not how to use it. See |.usage|_
for details on how to use |.FiPy|_.


.. list-table::
   :header-rows: 1
   
   * - 📝 Note
   * - It may be useful to set up a |.environment|_ before beginning
       the installation process.




-----------------------
Pre-Installed on Binder
-----------------------

A full |.FiPy|_ installation is available for basic exploration on
Binder_. The default notebook gives a rudimentary introduction to |.FiPy|_
syntax and, like any `Jupyter Notebook`_ interface, tab completion will help
you explore the package interactively.

.. _Binder:        https://mybinder.org/v2/gh/usnistgov/fipy/master
.. _Jupyter Notebook:    http://jupyter.org

.. _RECOMMENDED_METHOD:

------------------
Recommended Method
------------------

|CondaForge|_


.. list-table::
   :header-rows: 1
   
   * - 🔔️ Attention
   * - There are many ways to obtain the software
       packages necessary to run |.FiPy|_, but the most expedient way is
       with the conda_ package manager.  In addition to the scientific
       |.Python|_ stack, conda_ also provides virtual environment
       management.  Keeping separate installations is useful *e.g.* for
       comparing different versions of |.Python|_, or when
       the user does not have sufficient privileges to install software
       system-wide.

       In addition to the default packages, many other developers provide
       "channels" to distribute their own builds of a variety of software.
       These days, the most useful channel is conda-forge_, which provides
       everything necessary to install |.FiPy|_.


Install conda_
==============

`Install conda`_ or `install micromamba`_ on your computer.


.. _CREATE_CONDA_ENVIRONMENT:

Create a conda_ environment
===========================

Use one of the following methods to create a self-contained conda_
environment and then download and populate the environment with the
prerequisites for |.FiPy|_ from the conda-forge_ channel at
https://anaconda.org.  See `this discussion
<https://pythonspeed.com/articles/conda-dependency-management/>`_
of the merits of and relationship between the different methods.

* Conda_ environment files

  This option is the most upgradable in the future and probably the best
  for development.

  ::

    $ conda env create --name <MYFIPYENV> \
        --file environments/<SOLVER>-environment.yml


  .. list-table::
     :header-rows: 1
   
     * - 📝 Note
     * - You can try to include multiple solver suites using ``conda env
         update``, but be aware that different suites may have incompatible
         requirements, or may restrict installation to obsolete versions of
         Python.  Given that |.FiPy|_ can only use one solver suite during
         a run, installing more than one solver in an environment isn't
         necessary.


         .. list-table::
            :header-rows: 1
   
            * - 🔔️ Attention
            * - Successively updating an environment can be unpredictable, as later
                packages may conflict with earlier ones.  Unfortunately, ``conda
                env create`` `does not support multiple environment files
                <https://github.com/conda/conda/issues/9294>`_.

                Alternatively, combine the different
                :file:`environments/<SOLVER>-environment.yml` files you wish to
                use, along with `environment.yml` files for any other packages you
                are interested in (`conda-merge
                <https://github.com/amitbeka/conda-merge>`_ may prove useful).
                Then execute::

                  $ conda env create --name <MYFIPYENV> --file <MYMERGEDENVIRONMENT>.yml



* conda-lock_ lockfiles

  This option will pin all the packages, so is the most reproducible, but
  not particularly upgradable.  For most, this is the safest way to
  generate a FiPy environment that consistently works.

  ::

    $ conda-lock install --name <MYFIPYENV> \
        environments/locks/conda-<SOLVER>-lock.yml

  or, to be really explicit (and obviating the need for conda-lock_)::

    $ conda create --name <MYFIPYENV> \
        --file environments/locks/conda-<SOLVER>-<PLATFORM>.lock

* Directly from conda-forge_, picking and choosing desired packages

  This option is the most flexible, but has the highest risk of missing or
  incompatible packages.

  e.g.::

    $ conda create --name <MYFIPYENV> --channel conda-forge \
        python=3 numpy scipy matplotlib-base packaging mpich \
        mpi4py petsc4py mayavi "gmsh <4.0|>=4.5.2"

.. _conda-lock: https://github.com/conda/conda-lock




Install |.FiPy|_
====================

::

    $ conda install --name <MYFIPYENV> --channel conda-forge fipy


.. list-table::
   :header-rows: 1
   
   * - 📝 Note
   * - The `fipy conda-forge`_ package used to be "batteries included", but
       we found this to be too fragile.  It now only includes the bare
       minimum for |.FiPy|_ to function.


Enable conda_ environment
=========================

Enable your new environment with::

    $ conda activate <MYFIPYENV>

or::

    $ source activate <MYFIPYENV>

or, on Windows_::

    $ activate <MYFIPYENV>

You're now ready to move on to |.usage|_.


.. list-table::
   :header-rows: 1
   
   * - 📝 Note
   * - conda_ can be
       `quite <https://www.anaconda.com/blog/understanding-and-improving-condas-performance>`_
       `slow <https://medium.com/@marius.v.niekerk/conda-metachannel-f962241c9437>`_
       to resolve all dependencies when performing
       an installation.  You may wish to consider using the alternative
       mamba_ installation manager to speed things up.



.. list-table::
   :header-rows: 1
   
   * - 📝 Note
   * - On Linux_ and `Mac OS X`_, you should have a pretty complete system
       to run and visualize |.FiPy|_ simulations. On Windows_, there
       are fewer packages available via conda_, particularly amongst the
       sparse matrix |.solvers|_, but the system still should be
       functional.



.. list-table::
   :header-rows: 1
   
   * - 🔔️ Attention
   * - When installed via conda_ or |.pip|_, |.FiPy|_ will not include
       its |.part-colon-examples+examples|_.  These can be obtained by
       `cloning the repository`_ or downloading a `compressed archive`_.


.. _install conda: https://conda.io/projects/conda/en/latest/user-guide/install/
.. _install micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
.. _conda-forge: https://conda-forge.github.io/
.. _Mac OS X: http://www.apple.com/macosx/
.. _Linux: http://www.linux.org/
.. _Windows: http://www.microsoft.com/windows/
.. |CondaForge|    image:: https://anaconda.org/conda-forge/fipy/badges/version.svg
.. _CondaForge:    https://anaconda.org/conda-forge/fipy
.. _mamba: https://mamba.readthedocs.io/
.. _fipy conda-forge: https://anaconda.org/conda-forge/fipy


--------------
Obtaining FiPy
--------------

|.FiPy|_ is freely available for download via Git_ or as a
`compressed archive`_. Please see
|.documentation-colon-git|_ for instructions on obtaining |.FiPy|_
with Git_.


.. list-table::
   :header-rows: 1
   
   * - 🚩 Warning
   * - Keep in mind that if you choose to download the `compressed
       archive`_ you will then need to preserve your changes when upgrades
       to |.FiPy|_ become available (upgrades via Git_ will handle
       this issue automatically).


.. _Git:       https://github.com/usnistgov/fipy
.. _compressed archive:      https://github.com/usnistgov/fipy/releases

---------------
Installing FiPy
---------------

Details of the `Required Packages`_ and links are given below,
but for the courageous and the
impatient, |.FiPy|_ can be up and running quickly by simply
installing the following prerequisite packages on your system:

 * Python_

 * NumPy_

 * At least one of the |.solvers|_

 * At least one of the |.viewers|_ (|.FiPy|_'s tests will run
   without a viewer, but you'll want one for any practical work)

Other |.optionalpackages|_ add greatly to |.FiPy|_'s
capabilities, but are not necessary for an initial installation or to
simply run the test suite.

It is not necessary to formally install |.FiPy|_, but if you wish
to do so and you are confident that all of the requisite packages have
been installed properly, you can install it by typing::

    $ python -m pip install fipy

or by unpacking the archive and typing::

    $ python setup.py install

at the command line in the base |.FiPy|_ directory. You can also install
|.FiPy|_ in "development mode" by typing::

    $ python setup.py develop

which allows the source code to be altered in place and executed without
issuing further installation commands.

Alternatively, you may choose not to formally install |.FiPy|_ and
to simply work within the base directory instead. In this case or if you
are making a non-standard install (without admin privileges), read about
setting up your |.environment|_ before beginning the installation
process.

.. _REQUIREDPACKAGES:

-----------------
Required Packages
-----------------

.. warning:

   |.FiPy|_ will not run if the following items are not installed.

Python
======

http://www.python.org/

|.FiPy|_ is written in the |.Python|_ language and requires a
|.Python|_ installation to run. |.Python|_ comes pre-installed
on many operating systems, which you can check by opening a terminal
and typing ``python``, *e.g.*::

    $ python
    Python 3.10.10 | ...
    ...
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

If necessary, you can download_ and install it for your platform
<http://www.python.org/download>.


.. list-table::
   :header-rows: 1
   
   * - 📝 Note
   * - |.FiPy|_ `no longer supports Python 2 <sunset-python-2>`_.


.. _download: http://www.python.org/download/
.. _sunset-python-2: https://www.python.org/doc/sunset-python-2/

NumPy
=====

http://numpy.scipy.org

Obtain and install the |.NumPy|_ package. |.FiPy|_ requires at
least version 1.0 of NumPy_.

.. _OPTIONALPACKAGES:

-----------------
Optional Packages
-----------------

.. note:

    The following packages are not required to run |.FiPy|_, but they can
    be helpful.

Gmsh
====

http://www.geuz.org/gmsh/

|.Gmsh|_ is an application that allows the creation of irregular meshes.
When running in parallel, |.FiPy|_ requires a version of |.Gmsh|_
>= 2.5 and < 4.0 or >= 4.5.2.

SciPy
=====

http://www.scipy.org/

|.SciPy|_ provides a large collection of functions and tools that can
be useful for running and analyzing |.FiPy|_ simulations. Significantly
improved performance has been achieved with the judicious use of C language
inlining (see the |.flagsandenvironmentvariables|_ section for more
details), via the ``weave`` module.

.. note:

    A handful of test cases use functions from the |.SciPy|_
    library and will throw errors if it is missing.

------------------
Level Set Packages
------------------

To use the level set (:cite:`levelSetBook`) components of |.FiPy|_ one of the following is
required.

.. _SCIKITFMM:

Scikit-fmm
==========

http://packages.python.org/scikit-fmm/

Scikit-fmm_ is a python extension module which implements the fast
marching method.

.. _Scikit-fmm: http://packages.python.org/scikit-fmm/

.. _LSMLIBDOC:

LSMLIB
======

http://ktchu.serendipityresearch.org/software/lsmlib/index.html

The Level Set Method Library (LSMLIB_) provides support for the serial
and parallel simulation of implicit surface and curve dynamics in two-
and three-dimensions.

Install LSMLIB_ as per the instructions on the website. Additionally
PyLSMLIB_ is required. To install, follow the instructions on the
website,
https://github.com/ktchu/LSMLIB/tree/master/pylsmlib#pylsmlib.

.. _PyLSMLIB: https://github.com/ktchu/LSMLIB/tree/master/pylsmlib#pylsmlib
.. _LSMLIB: http://ktchu.serendipityresearch.org/software/lsmlib/index.html

.. _ENVIRONMENT:

-----------------------
Development Environment
-----------------------

It is often preferable to not formally install packages in the system
directories. The reasons for this include:

 * developing or altering the package source code,

 * trying out a new package along with its dependencies without
   violating a working system,

 * dealing with conflicting packages and dependencies,

 * or not having admin privileges.

To avoid tampering with the system Python_ installation, you can employ one
of the utilities that manage packages and their dependencies independently
of the system package manager and the system directories.  These utilities
include conda_, Nix_, Stow_, Virtualenv_ and Buildout_, amongst others.
Conda_ and Nix_ are only ones of these we have the resources to support.

|.create_conda_environment|_ for development, followed by::

   $ source activate <MYFIPYENV>
   $ python -m pip install scikit-fmm
   $ git clone https://github.com/usnistgov/fipy.git
   $ cd fipy
   $ python setup.py develop

.. _Conda: https://conda.io
.. _Stow: http://savannah.gnu.org/projects/stow/
.. _Buildout: http://pypi.python.org/pypi/zc.buildout
.. _Virtualenv: https://virtualenv.pypa.io

.. _documentation:GIT:

---------
Git usage
---------

All stages of |.FiPy|_ development are archived in a Git
repository at GitHub_. You can browse through the code at
https://github.com/usnistgov/fipy and, using a `Git client`_, you can
download various tagged revisions of |.FiPy|_ depending on your needs.


.. list-table::
   :header-rows: 1
   
   * - 🔔️ Attention
   * - Be sure to follow |.installation|_ to obtain all the prerequisites for
       |.FiPy|_.


Git client
==========

A ``git`` client application is needed in order to fetch files from our
repository. This is provided on many operating systems (try executing
``which git``) but needs to be installed on many others. The sources to
build Git, as well as links to various pre-built binaries for
different platforms, can be obtained from http://git-scm.com/.

Git branches
============

In general, most users will not want to download the very latest state of
|.FiPy|_, as these files are subject to active development and may not behave
as desired. Most users will not be interested in particular version numbers
either, but instead with the degree of code stability. Different branches are
used to indicate different stages of |.FiPy|_ development. For the
most part, we follow `a successful Git branching model`_. You will
need to decide on your own risk tolerance when deciding which stage of
development to track.

.. _cloning the repository:

A fresh copy of the |.FiPy|_ source code  can be obtained with::

   $ git clone https://github.com/usnistgov/fipy.git

An existing Git checkout of FiPy can be shifted to a different `<branch>` of
development by issuing the command::

   $ git checkout <branch>

in the base directory of the working copy. The main branches for FiPy are:

``master``
    designates the (ready to) release state of FiPy. This code is stable
    and should pass all of the tests (or should be documented that it does
    not).

Past releases of FiPy are tagged as

``x.y.z``
    Any released version of FiPy will be designated with a fixed tag: The
    current version of FiPy is |version|.  (Legacy ``version-x_y_z`` tags
    are retained for historical purposes, but won't be added to.)

Tagged releases can be found with::

   $ git tag --list

Any other branches will not generally be of interest to most users.


.. list-table::
   :header-rows: 1
   
   * - 📝 Note
   * - For some time now, we have done all significant development work on
       branches, only merged back to ``master`` when the tests pass
       successfully.  Although we cannot guarantee that ``master`` will never
       be broken, you can always check our |.continuousintegration|_ status
       to find the most recent revision that it is running acceptably.

       Historically, we merged to ``develop`` before merging to ``master``.  We
       no longer do this, although for time being, ``develop`` is kept
       synchronized with ``master``.  In a future release, we will remove the
       ``develop`` branch altogether.


For those who are interested in learning more about Git, a wide variety of
online sources are available, starting with the `official Git website`_.
The `Pro Git book`_ :cite:`ProGit` is particularly instructive.

.. _official Git website: http://git-scm.com/

.. _Pro Git book: http://git-scm.com/book

.. _GitHub: https://github.com/usnistgov/fipy

.. _a successful Git branching model: http://nvie.com/posts/a-successful-git-branching-model/


---
Nix
---

.. _nixinstall:

Nix Installation
================

|.FiPy|_ now has a `Nix`_ expression for installing |.FiPy|_
using `Nix`_. `Nix`_ is a powerful package manager for Linux and other
Unix systems that makes package management reliable and
reproducible. The recipe works on both Linux and Mac OS X. Go to
`nix.dev`_ to get started with Nix.

Installing
----------

Once you have a working Nix installation use::

    $ nix develop

in the base |.FiPy|_ directory to install |.FiPy|_ with Python
3 by default. ``nix develop`` drops the user into a shell with a working
version of |.FiPy|_. To test your installation use::

    $ nix develop --command bash -c "python setup.py test"


.. list-table::
   :header-rows: 1
   
   * - 📝 Note
   * - The SciPy solvers are the only available solvers currently.



.. _Nix: https://nixos.org/nix/
.. _Nixpkgs:  https://nixos.org/nixpkgs/
.. _nix.dev: https://nix.dev

