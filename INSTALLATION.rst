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
      this page to fipy@nist.gov.

-----------------------
Pre-Installed on Binder
-----------------------

A full :term:`FiPy` installation is available for basic exploration on
Binder_. The default notebook gives a rudimentary introduction to :term:`FiPy`
syntax and, like any `Jupyter Notebook`_ interface, tab completion will help
you explore the package interactively.

.. _Binder:        https://mybinder.org/v2/gh/usnistgov/fipy/develop
.. _Jupyter Notebook:    http://jupyter.org

------------------
Recommended Method
------------------

.. attention::

   There are many ways (described further down) to obtain the software
   packages necessary to run :term:`FiPy`, but the most expedient way is
   with the conda_ package manager.

  * `install Miniconda`_ on your computer 
  * run::

      $ conda create --name <MYFIPYENV> --channel guyer --channel conda-forge fipy nomkl

    .. note::

       This command creates a self-contained conda_ environment and then
       downloads and populates the environment with the prerequisites for
       :term:`FiPy` from the channels guyer_ and conda-forge_ at
       https://anaconda.org.

    .. attention::

       Note, this does not work on Windows. On that platform, you should be able to do::

        conda create --name <MYFIPYENV> --channel guyer --channel conda-forge python numpy scipy matplotlib mayavi
        activate <MYFIPYENV>
        pip install fipy

       There are presently no conda packages of any solver suite but scipy available for Windows.


  * enable this new environment with::

      $ source activate <MYFIPYENV>

    .. note::

       ``$ activate <MYFIPYENV>`` on Windows_

  * move on to :ref:`USAGE`. 

    .. note::

       On Linux_ and `Mac OS X`_, you should have a pretty complete system
       to run and visualize :term:`FiPy` simulations. On Windows_, there
       are fewer packages available via conda_, particularly amongst the
       sparse matrix :ref:`SOLVERS`, but the system still should be
       functional.

.. _install Miniconda: http://conda.pydata.org/docs/install/quick.html
.. _guyer: https://anaconda.org/guyer
.. _conda-forge: https://conda-forge.github.io/

--------------------------
Installing Python Packages
--------------------------

In general, it is best to use the following order of precedence when
installing packages:

 * Use the operating system package manager, if possible.
 * Use the conda_ package management system, which handles both 
   :term:`Python` and non-:term:`Python` packages and provides facilities
   for self-contained environments with different combinations of
   :term:`Python` packages, libraries, and applications.
 * Use the 
   `pip installs python <http://www.pip-installer.org/>`_ 
   (:term:`pip`) tool to obtain software from the 
   `Python Package Index <http://pypi.python.org/pypi>`_ 
   (:term:`PyPI`) repository::

     $ pip install package

   .. warning::

      :term:`pip` takes care of dependencies that are themselves
      :term:`Python` packages. It does not deal with non-:term:`Python`
      dependencies.

 * Download the packages manually, unpack and run::

     $ python setup.py install

   Further information about each ``setup.py`` script is available by
   typing::

     $ python setup.py --help

Many of the packages listed below have prebuilt installers for
different platforms (particularly for Windows).  These installers can
save considerable time and effort compared to configuring and building
from source, although they frequently comprise somewhat older versions
of the respective code.  Whether building from source or using a
prebuilt installer, please read and follow explicitly any instructions
given in the respective packages' :file:`README` and
:file:`INSTALLATION` files.

--------------
Obtaining FiPy
--------------

:term:`FiPy` is freely available for download via Git_ or as a
compressed archive from
<http://www.ctcms.nist.gov/fipy/download>. Please see
:ref:`documentation:GIT` for instructions on obtaining :term:`FiPy`
with Git_.

.. warning::

   Keep in mind that if you choose to download the `compressed
   archive`_ you will then need to preserve your changes when upgrades
   to :term:`FiPy` become available (upgrades via Git_ will handle
   this issue automatically).

.. _Git:       https://github.com/usnistgov/fipy
.. _compressed archive:      http://www.ctcms.nist.gov/fipy/download/

---------------
Installing FiPy
---------------

Details of the `Required Packages`_ and links are given below and in
`platform-specific instructions`_, but for the courageous and the
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

    $ pip install fipy

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
    Python 2.3 (#1, Sep 13 2003, 00:49:11) 
    ...
    Type "help", "copyright", "credits" or "license" for more information.
    >>> 

If necessary, you can download_ and install it for your platform
<http://www.python.org/download>.

.. note::

   :term:`FiPy` requires at least version 2.4.x of :term:`Python`. See
   the specialized instructions if you plan on :ref:`RunningUnderPython3`.

.. _download: http://www.python.org/download/

:term:`Python` along with many of :term:`FiPy`'s required and optional
packages is available with one of the following distributions.

conda
-----

http://conda.pydata.org

This package manager provides a wide array of both :term:`Python`-based and
general scientific packages. In addition to the default packages, many
other developers (including us) provide "channels" to distribute their own
builds of a variety of software.

In a given conda_ environment, you can install :term:`FiPy` with::

    $ conda install --channel guyer --channel conda-forge fipy


.. _EPD:

Scientific Python Packages
==========================

:term:`FiPy` depends on or benefits from the so-called *scientific Python*
stack, which includes :term:`IPython`, :term:`Matplotlib`, :term:`NumPy`, 
:term:`pandas`, and :term:`SciPy`. These packages provide Python with advanced
numerical and graphical capabilities, important for the analysis and visual
representation of scientific data. The following :term:`Python` distributions
provide the :term:`Python` interpreter, scientific Python stacks, and searchable
package management for all three major operating systems:
Linux_, `Mac OS X`_, and Windows_.

.. attention::

    :term:`Python` distributions include cryptographic software packages,
    and may be subject to export control laws.

.. _PYPI:

Pip Installs Python
-------------------

If you have the reference :term:`Python` distribution installed, :term:`pip`
provides a straightforward way to install and manage the scientific :term:`Python`
stack. :term:`pip` interfaces with :term:`PyPI`, searchable on the Web or through
the command line by opening a terminal and typing::

    $ pip search <package name>

:term:`PyPI` packages are maintained by the individual authors, so installations
using :term:`pip` run a somewhat elevated risk of unmet dependencies as the
independent codes evolve. Stable releases of :term:`FiPy` are packaged for
distribution through :term:`PyPI`.

.. _ANACONDA:

Continuum Analytics Anaconda
----------------------------

http://continuum.io/anaconda

In addition to the scientific :term:`Python` stack, the Anaconda package manager
also provides virtual environment management. Keeping separate installations is useful
*e.g.* for comparing :term:`Python` 2 and :term:`Python` 3 software stacks, or
when the user does not have sufficient provileges to install software system-wide.

.. attention::

   :term:`PySparse` and :term:`FiPy` are not presently included in Anaconda,
   so you will need to separately install them manually.

.. _ECP:

Enthought Canopy
----------------

http://www.enthought.com/products/canopy

In addition to the core scientific :term:`Python` stack, Canopy includes packages
for a very large number of software projects. Canopy supports the venv virtual
environment manager.

.. attention::

   :term:`PySparse` and :term:`FiPy` are not presently included in Canopy,
   so you will need to separately install them manually.

.. _PYTHONXY:

Python(x,y)
-----------

http://python-xy.github.io

A non-commercial scientific :term:`Python` distribution developed specifically
for Windows_.

.. attention::

   :term:`PySparse` and :term:`FiPy` are not presently included in
   Python(x,y), so you will need to separately install them manually.

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

------------------------------
Platform-Specific Instructions
------------------------------

:term:`FiPy` is `tested regularly`_ on `Mac OS X`_, `Debian Linux`_,
`Ubuntu Linux`_, and `Windows XP`_. We welcome reports of compatibility
with other systems, particularly if any additional steps are necessary to
install (see `Miscellaneous Build Recipes`_ for user contributed
installation tips).

The only elements of :term:`FiPy` that are likely to be
platform-dependent are the :ref:`VIEWERS` but at least one viewer
should work on each platform.  All other aspects should function on
any platform that has a recent :term:`Python` installation.

Mac OS X Installation
=====================

There is no official package manager for `Mac OS X`_, but there are 
several third-party package managers that provide many, but not all of 
:term:`FiPy`'s :ref:`REQUIREDPACKAGES` and :ref:`OPTIONALPACKAGES`. 
Options include

 Fink_
   is based on the Debian package management system. It installs all of its
   dependencies into :file:`/sw`.
 MacPorts_
   is a package manager originally part of OpenDarwin. It installs all of 
   its dependencies into :file:`/opt`.
 Homebrew_
   is a recent, lightweight package manager based on Ruby scripts. It
   installs all of its dependencies into :file:`/usr/local` (although it
   can be directed not to).

In addition, there is an :ref:`ECP` installer for `Mac OS X`_.

.. attention::

   :term:`PySparse` and :term:`FiPy` are not presently included in any of
   these package managers or installers, so you will need to separately
   install them manually.

We presently find that the combination of Homebrew_ and :term:`pip` is a 
pretty straightforward way to get most of :term:`FiPy`'s prerequesites. 
See the `Miscellaneous Build Recipes`_ for up-to-date directions.

.. _Fink: http://www.finkproject.org/
.. _MacPorts: http://www.macports.org/
.. _Homebrew: http://mxcl.github.com/homebrew/

Windows Installation
====================

There is no official package manager for Windows_, but the :ref:`ECP` 
and :ref:`PYTHONXY` installers provide most of :term:`FiPy`'s 
prerequisites.

.. attention::

   :term:`PySparse` and :term:`FiPy` are not presently included in Canopy or
   Python(x,y), so you will need to separately install them manually.

Ubuntu/Debian Installation
==========================

:term:`FiPy` now has a `.deb` for Ubuntu/Debian systems that can be
downloaded from <http://www.ctcms.nist.gov/fipy/download>. Simply run::

  $ VERSION=x.y-z  # choose the version you want
  $ apt-get install gmsh libsuperlu3 python-central python-sparse
  $ curl -O http://www.ctcms.nist.gov/fipy/download/python-fipy_${VERSION}_all.deb
  $ dpkg -i python-fipy_${VERSION}_all.deb

to install. The `.deb` includes dependencies for all of the
:ref:`REQUIREDPACKAGES` and :ref:`OPTIONALPACKAGES`.

.. _tested regularly: http://matforge.org/fipy/build
.. _Mac OS X: http://www.apple.com/macosx/
.. _Linux: http://www.linux.org/
.. _Debian Linux: http://www.debian.org/
.. _RedHat Linux: http://www.redhat.com/
.. _OpenSUSE Linux: http://www.opensuse.org/
.. _Ubuntu Linux: http://www.ubuntu.com/
.. _Solaris: http://oracle.com/solaris
.. _Windows: http://www.microsoft.com/windows/
.. _Windows XP: http://www.microsoft.com/windowsxp/

Miscellaneous Build Recipes
===========================

We often post miscellaneous installation instructions on the
:term:`FiPy` blog_ and wiki_ pages. The most useful of these include:

 * `Installing FiPy on Mac OS X using Homebrew`_
 * `Building a 64-bit scientific python environment for FiPy from source`_
 * `Installing FiPy with pip`_

.. note::

    We encourange you to contribute your own build recipes on the wiki_
    if they are significantly different.

.. _Installing FiPy on Mac OS X using Homebrew: http://matforge.org/fipy/wiki/InstallFiPy/MacOSX/HomeBrew
.. _Building a 64-bit scientific python environment for FiPy from source: http://matforge.org/fipy/wiki/InstallFiPy/MacOSX/SnowLeopard
.. _Installing FiPy with pip: http://matforge.org/fipy/wiki/InstallFiPy/PipInstallsPython
.. _wiki: http://matforge.org/fipy
.. _blog: http://matforge.org/fipy/blog

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

The simplest way to use a :term:`Python` package without installing it
is to work in the base directory of the unpacked package and set 
the :envvar:`PYTHONPATH` environment variable to "``.``". In order to work in an
directory other than the package's base directory, the :envvar:`PYTHONPATH`
environment variable must be set to ":file:`~/path/to/package`". This
method of working is adequate for one package, but quickly becomes
unmanageable with multiple :term:`Python` packages. In order to manage
multiple packages, it is better to choose a standard location other
than the default installation path. 

If you do not have administrative privileges on your computer, or if
for any reason you don't want to tamper with your existing
:term:`Python` installation, most packages (including :term:`FiPy`)
will allow you to install to an alternative location.  Instead of
installing these packages with ``python setup.py install``, you would
use :samp:`python setup.py install --home={dir}`, where :samp:`{dir}`
is the desired installation directory (usually "``~``" to indicate
your home directory).  You will then need to append
:file:`{dir}/lib/python` to your :envvar:`PYTHONPATH` environment
variable.  See the `Alternate Installation`_ section of the
:term:`Python` document "`Installing Python Modules`_"
:cite:`InstallingPythonModules` for more information, such as circumstances
in which you should use :option:`--prefix` instead of
:option:`--home`.

.. _Alternate Installation: http://docs.python.org/inst/alt-install-windows.html

.. _Installing Python Modules: http://docs.python.org/inst/

An alternative to setting the :envvar:`PYTHONPATH` is to employ one of the
utilities that manage packages and their dependencies independently of
the system package manager and the system directories. These utilities
include Conda_, Stow_, Virtualenv_ and zc.buildout_, amongst others. Here we'll
describe the use of Conda_, which we highly recommend.

.. _Stow: http://savannah.gnu.org/projects/stow/
.. _zc.buildout: http://pypi.python.org/pypi/zc.buildout
.. _Virtualenv: https://virtualenv.pypa.io
.. _Conda: https://conda.io

.. _documentation:GIT:

.. include:: documentation/GIT.rst
