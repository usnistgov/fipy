========
Overview
========

.. only:: latex

   :term:`FiPy` is an object oriented, partial differential equation (PDE) solver,
   written in Python_ [Python]_, based on a standard finite volume (FV)
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

   The :term:`FiPy` framework includes terms for transient diffusion,
   convection and standard sources, enabling the solution of arbitrary
   combinations of coupled elliptic, hyperbolic and parabolic PDEs. Currently
   implemented models include phase field
   [BoettingerReview:2002]_ [ChenReview:2002]_ [McFaddenReview:2002]_ treatments of
   polycrystalline, dendritic, and electrochemical phase transformations as
   well as a level set treatment of the electrodeposition process
   [NIST:damascene:2001]_.

.. only:: latex
  
   The latest information about :term:`FiPy` can be found at
   http://www.ctcms.nist.gov/fipy/.

---------------------------------
Even if you don't read manuals...
---------------------------------

...please read :ref:`INSTALLATION` and the :ref:`FAQ`. 

--------------------------------
What's new in version |VERSION|?
--------------------------------

.. warning::

   :term:`FiPy` 2 brings unavoidable syntax changes. Please see
   :mod:`examples.updating.update1_0to2_0` for guidance on the changes that
   you will need to make to your :term:`FiPy` 1.x scripts.

The significant changes since version 1.2 are:

- :class:`CellVariable` and :class:`FaceVariable` objects can hold values of any 
  rank.

- Much simpler syntax for specifying :class:`Cells`\s for initial conditions and 
  :class:`Face`\s for boundary conditions.

- Automated determination of the Peclet number and partitioning of 
  :class:`ImplicitSourceTerm` coefficients between the matrix diagonal and the
  right-hand-side-vector.

- Simplified :class:`Viewer` syntax.

- Support for the `Trilinos solvers`_.

- Support for anisotropic diffusion coefficients.

.. _Trilinos solvers: http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://trilinos.sandia.gov

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

You can communicate with the :term:`FiPy` developers and with other users via our
`mailing list`_ and we welcome you to use the `tracking
system`_ for bugs, support requests, feature requests and
patch submissions <http://matforge.org/fipy/report>. We welcome collaborative efforts on this project.

:term:`FiPy` is a member of MatForge_, a project of the `Materials Digital Library
Pathway`_. This National Science Foundation funded service provides
management of our public source code repository, our bug tracking system, and
a "wiki" space for public contributions of code snippets, discussions, and
tutorials.

.. toctree::

   documentation/MAIL

------------------------
Conventions and Notation
------------------------

:term:`FiPy` is driven by Python_ script files than you can view or modify in any
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
current working directory is the :term:`FiPy` distribution directory, refered to
as the "base directory", such that::

    examples/diffusion/steadyState/mesh1D.py

will correspond to, *e.g.*::

    /some/where/FiPy-X.Y/examples/diffusion/steadyState/mesh1D.py

Paths will always be rendered using POSIX conventions (path elements
separated by "``/``").  Any references of the form::

    examples.diffusion.steadyState.mesh1D

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

.. include:: documentation/VERSION.txt

