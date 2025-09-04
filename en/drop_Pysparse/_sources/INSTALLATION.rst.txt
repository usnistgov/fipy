.. _INSTALLATION:

============
Installation
============

The :term:`FiPy` finite volume PDE solver relies on several
third-party packages.  It is *best to obtain and install those first*
before attempting to install :term:`FiPy`. This document explains how
to install :term:`FiPy`, not how to use it. See :ref:`USAGE`
for details on how to use :term:`FiPy`.

.. note::

   It may be useful to set up a :ref:`ENVIRONMENT` before beginning
   the installation process.

.. only:: html

   .. note::

      By selecting the links on this page, you will be leaving NIST
      webspace. We have provided these links to other web sites because
      they may have information that would be of interest to you. No
      inferences should be drawn on account of other sites being
      referenced, or not, from this page. There may be other web sites that
      are more appropriate for your purpose. NIST does not necessarily
      endorse the views expressed, or concur with the facts presented on
      these sites. Further, NIST does not endorse any commercial products
      that may be mentioned on these sites. Please address comments about
      this page to fipy@list.nist.gov.

-----------------------
Pre-Installed on Binder
-----------------------

A full :term:`FiPy` installation is available for basic exploration on
Binder_. The default notebook gives a rudimentary introduction to :term:`FiPy`
syntax and, like any `Jupyter Notebook`_ interface, tab completion will help
you explore the package interactively.

.. _Binder:        https://mybinder.org/v2/gh/usnistgov/fipy/master
.. _Jupyter Notebook:    http://jupyter.org

.. _RECOMMENDED_METHOD:

------------------
Recommended Method
------------------

|CondaForge|_

.. attention::

   There are many ways to obtain the software
   packages necessary to run :term:`FiPy`, but the most expedient way is
   with the conda_ package manager.  In addition to the scientific
   :term:`Python` stack, conda_ also provides virtual environment
   management.  Keeping separate installations is useful *e.g.* for
   comparing :term:`Python` 2 and :term:`Python` 3 software stacks, or when
   the user does not have sufficient privileges to install software
   system-wide.

   In addition to the default packages, many other developers provide
   "channels" to distribute their own builds of a variety of software.
   These days, the most useful channel is conda-forge_, which provides
   everything necessary to install :term:`FiPy`.

Install conda_
==============

`Install conda`_ or `install micromamba`_ on your computer.


.. include:: ../../environments/README.rst


Install :term:`FiPy`
====================

::

    $ conda install --name <MYFIPYENV> --channel conda-forge fipy

.. note::

   The `fipy conda-forge`_ package used to be "batteries included", but
   we found this to be too fragile.  It now only includes the bare
   minimum for :term:`FiPy` to function.

Enable conda_ environment
=========================

Enable your new environment with::

    $ conda activate <MYFIPYENV>

or::

    $ source activate <MYFIPYENV>

or, on Windows_::

    $ activate <MYFIPYENV>

You're now ready to move on to :ref:`USAGE`.

.. note::

   conda_ can be
   `quite <https://www.anaconda.com/blog/understanding-and-improving-condas-performance>`_
   `slow <https://medium.com/@marius.v.niekerk/conda-metachannel-f962241c9437>`_
   to resolve all dependencies when performing
   an installation.  You may wish to consider using the alternative
   mamba_ installation manager to speed things up.

.. note::

   On Linux_ and `Mac OS X`_, you should have a pretty complete system
   to run and visualize :term:`FiPy` simulations. On Windows_, there
   are fewer packages available via conda_, particularly amongst the
   sparse matrix :ref:`SOLVERS`, but the system still should be
   functional. Significantly, you will need to download and install
   :term:`Gmsh` manually when using Python 2.7.

.. attention::

   When installed via conda_ or :term:`pip`, :term:`FiPy` will not include
   its :ref:`examples <part:examples>`.  These can be obtained by
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

:term:`FiPy` is freely available for download via Git_ or as a
`compressed archive`_. Please see
:ref:`documentation:GIT` for instructions on obtaining :term:`FiPy`
with Git_.

.. warning::

   Keep in mind that if you choose to download the `compressed
   archive`_ you will then need to preserve your changes when upgrades
   to :term:`FiPy` become available (upgrades via Git_ will handle
   this issue automatically).

.. _Git:       https://github.com/usnistgov/fipy
.. _compressed archive:      https://github.com/usnistgov/fipy/releases

---------------
Installing FiPy
---------------

Details of the `Required Packages`_ and links are given below,
but for the courageous and the
impatient, :term:`FiPy` can be up and running quickly by simply
installing the following prerequisite packages on your system:

 * Python_

 * NumPy_

 * At least one of the :ref:`SOLVERS`

 * At least one of the :ref:`VIEWERS` (:term:`FiPy`'s tests will run
   without a viewer, but you'll want one for any practical work)

Other :ref:`OPTIONALPACKAGES` add greatly to :term:`FiPy`'s
capabilities, but are not necessary for an initial installation or to
simply run the test suite.

It is not necessary to formally install :term:`FiPy`, but if you wish
to do so and you are confident that all of the requisite packages have
been installed properly, you can install it by typing::

    $ python -m pip install fipy

or by unpacking the archive and typing::

    $ python setup.py install

at the command line in the base :term:`FiPy` directory. You can also install
:term:`FiPy` in "development mode" by typing::

    $ python setup.py develop

which allows the source code to be altered in place and executed without
issuing further installation commands.

Alternatively, you may choose not to formally install :term:`FiPy` and
to simply work within the base directory instead. In this case or if you
are making a non-standard install (without admin privileges), read about
setting up your :ref:`ENVIRONMENT` before beginning the installation
process.

.. _REQUIREDPACKAGES:

-----------------
Required Packages
-----------------

.. warning:

   :term:`FiPy` will not run if the following items are not installed.

Python
======

http://www.python.org/

:term:`FiPy` is written in the :term:`Python` language and requires a
:term:`Python` installation to run. :term:`Python` comes pre-installed
on many operating systems, which you can check by opening a terminal
and typing ``python``, *e.g.*::

    $ python
    Python 2.7.15 | ...
    ...
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

If necessary, you can download_ and install it for your platform
<http://www.python.org/download>.

.. note::

   :term:`FiPy` requires at least version 2.7.x of :term:`Python`.

.. _download: http://www.python.org/download/

:term:`Python` along with many of :term:`FiPy`'s required and optional
packages is available with one of the following distributions.

NumPy
=====

http://numpy.scipy.org

Obtain and install the :term:`NumPy` package. :term:`FiPy` requires at
least version 1.0 of NumPy_.

.. _OPTIONALPACKAGES:

-----------------
Optional Packages
-----------------

.. note:

    The following packages are not required to run :term:`FiPy`, but they can
    be helpful.

Gmsh
====

http://www.geuz.org/gmsh/

:term:`Gmsh` is an application that allows the creation of irregular meshes.
When running in parallel, :term:`FiPy` requires a version of :term:`Gmsh`
>= 2.5 and < 4.0 or >= 4.5.2.

SciPy
=====

http://www.scipy.org/

:term:`SciPy` provides a large collection of functions and tools that can
be useful for running and analyzing :term:`FiPy` simulations. Significantly
improved performance has been achieved with the judicious use of C language
inlining (see the :ref:`FlagsAndEnvironmentVariables` section for more
details), via the :mod:`weave` module.

.. note:

    A handful of test cases use functions from the :term:`SciPy`
    library and will throw errors if it is missing.

------------------
Level Set Packages
------------------

To use the level set (:cite:`levelSetBook`) components of :term:`FiPy` one of the following is
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

:ref:`CREATE_CONDA_ENVIRONMENT` for development, followed by::

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

.. include:: GIT.rst

---
Nix
---

.. _nixinstall:

.. include:: NIX-README.rst
