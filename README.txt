========
Overview
========

|FiPy| is an object oriented, partial differential equation (PDE) solver,
written in Python_ |citePython|, based on a standard finite volume (FV)
approach.  The framework has been developed in the `Metallurgy Division`_
and Center for Theoretical and Computational Materials Science (CTCMS_), in
the Materials Science and Engineering Laboratory (MSEL_) at the National
Institute of Standards and Technology (NIST_).

The solution of coupled sets of PDEs is ubiquitous to the numerical
simulation of science problems.  Numerous PDE solvers exist, using a
variety of languages and numerical approaches. Many are proprietary,
expensive and difficult to customize.  As a result, scientists spend
considerable resources repeatedly developing limited tools for
specific problems.  Our approach, combining the FV method and Python_,
provides a tool that is extensible, powerful and freely available. A
significant advantage to Python_ is the existing suite of tools for
array calculations, sparse matrices and data rendering.

The |FiPy| framework includes terms for transient diffusion, convection and
standard sources, enabling the solution of arbitrary combinations of
coupled elliptic, hyperbolic and parabolic PDEs.  Currently implemented
models include phase field |citePhaseField| treatments of polycrystalline,
dendritic, and electrochemical phase transformations as well as a level set
treatment of the electrodeposition process |citeCEAC|.

.. raw:: html
  
   <!--

The latest information about |FiPy| can be found at
http://www.ctcms.nist.gov/fipy/.

.. raw:: html
  
   -->

---------------------------------
Even if you don't read manuals...
---------------------------------

...please read |INSTALLATION-txt| and |the FAQ|. 

-----------------------------------
What's new in version |VERSION|
-----------------------------------

The significant changes since version 1.0 are:

- Memory efficiency has been improved in a number of ways, but most
  significantly by:

  * not caching all intermediate ``Variable`` values.
  * introducing ``UniformGrid`` classes that calculate geometric
    arrays on the fly.

  Details of these improvements are presented in |EFFICIENCY-txt|.

- Installation on Windows has been made considerably easier by
  constructing executable installers for |FiPy| and its
  dependencies. Instructions for Windows installation can be found in
  |WINDOWS-INSTALLATION-txt|.

- The arithmetic for ``Variable`` subclasses now works, and returns
  sensible answers. For example, ``VectorCellVariable * CellVariable``
  returns a ``VectorCellVariable``.

- ``PeriodicGrid`` meshes have been implemented. Currently, however,
  there and no examples of their use in the manual.

- Many of the examples have been completely rewritten

  * A basic 1D diffusion problem now serves as a general tutorial for 
    setting up any problem in |FiPy|. 
  * Several more phase field examples have been added that should make it 
    clearer how to get from the simple 1D case to the more elaborate 
    multicomponent, multidimensional, and anisotropic models.
  * The "Superfill" examples have been substantially improved with better
    functionality and documentation.
  * An example of fluid flow with the classic Stokes moving lid has been 
    added.

- A clear distinction has been made between solving an equation via `solve()`
  and iterating an non-linear equation to solution via `sweep()`. An extensive 
  explanation of the concepts involved has been added to |the FAQ|.

- Added a `MultiViewer` class that automatically groups several viewers 
  together if the variables couldn't be displayed by a single viewer.

- The abbreviated syntax ``from fipy import Class`` or ``from fipy import *``
  promised in version 1.0 actually works now. The examples all still use the
  fully qualified names.

- The repository has been converted from a CVS to a Subversion_
  repository. Details on how to check out the new repository are given
  in the |INSTALLATION-txt|.

- The |FiPy| repository has also been moved from Sourceforge_ to 
  MatForge_.

-------------------------
Download and Installation
-------------------------

Please refer to |INSTALLATION-txt| for details on download and
installation. |FiPy| can be redistributed and/or modified
freely, provided that any derivative works bear some notice that they
are derived from it, and any modified versions bear some notice that
they have been modified.

-------
Support
-------

You can communicate with the |FiPy| developers and with other users via our
`mailing list`_ and we welcome you to use the `tracking
system`_ for bugs, support requests, feature requests and
patch submissions |citeMailingTracking|. We welcome collaborative efforts on this project.

|FiPy| is a member of MatForge_, a project of the `Materials Digital Library
Pathway`_. This National Science Foundation funded service provides
management of our public source code repository, our bug tracking system, and
a "wiki" space for public contributions of code snippets, discussions, and
tutorials.

------------------------
Conventions and Notation
------------------------

|FiPy| is driven by Python_ script files than you can view or modify in any
text editor.  |FiPy| sessions are invoked from a command-line shell, such
as ``tcsh`` or ``bash``.

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

are intended to indicate an interactive session in the Python_ interpreter.
We will refer to these as "interactive sessions" or as "doctest blocks".
The text "``>>>``" at the beginning of a line denotes the *primary prompt*,
calling for input of a Python_ command.  The text "``...``" denotes the
*secondary prompt*, which calls for input that continues from the line
above, when required by Python_ syntax.  All remaining lines, which begin
at the left margin, denote output from the Python_ interpreter.  In all
cases, the prompt is supplied by the Python_ interpreter and should not be
typed by you.

.. warning::

   Python_ is sensitive to indentation and care should be taken to enter
   text exactly as it appears in the examples.

When references are made to file system paths, it is assumed that the
current working directory is the |FiPy| distribution directory, refered to
as the "base directory", such that::

    examples/diffusion/steadyState/mesh1D/input.py

will correspond to, *e.g.*::

    /some/where/FiPy-1.0/examples/diffusion/steadyState/mesh1D/input.py

Paths will always be rendered using POSIX conventions (path elements
separated by "``/``").  Any references of the form::

    examples.diffusion.steadyState.mesh1D.input

are in the Python_ module notation and correspond to the equivalent POSIX
path given above.

We may at times use a 

.. note::

   to indicate something that may be of interest

or a

.. warning::

   to indicate something that could cause serious problems.

.. _MSEL:                 http://www.msel.nist.gov/
.. _CTCMS:                http://www.ctcms.nist.gov/
.. _Metallurgy Division:  http://www.metallurgy.nist.gov/
.. _NIST:                 http://www.nist.gov/
.. _Python:               http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://www.python.org/
.. _Subversion:           http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://matforge.org/fipy/browser
.. _compressed archive:   http://www.ctcms.nist.gov/fipy/download/FiPy-1.1.tar.gz
.. _tracking system:      http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://matforge.org/fipy/report
.. _mailing list:         http://www.ctcms.nist.gov/fipy/mail.html
.. _Sourceforge:          http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://www.sourceforge.net/projects/fipy
.. _Materials Digital Library Pathway: http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://matdl.org
.. _MatForge:             http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://matforge.org/

.. include:: utils/include.txt
.. include:: documentation/VERSION.txt

.. |FiPy| replace:: |htmlFiPy| |latexFiPy|
.. |INSTALLATION-txt| replace:: |htmlINSTALL| |latexINSTALL|
.. |EFFICIENCY-txt| replace:: |htmlEFFICIENCY| |latexEFFICIENCY|
.. |the FAQ| replace:: |htmlFAQ| |latexFAQ|
.. |WINDOWS-INSTALLATION-txt| replace:: |htmlWINDOWS-INSTALLATION| |latexWINDOWS-INSTALLATION|

.. |citePython| raw:: latex

   \cite{Python}

.. |citePhaseField| raw:: latex

   \cite{BoettingerReview:2002,ChenReview:2002,McFaddenReview:2002}

.. |citeCEAC| raw:: latex

   \cite{NIST:damascene:2001}

.. |citeMailingTracking| raw:: latex

   \cite{FiPyMailingList,FiPyBugTracker}

