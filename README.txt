========
Overview
========

|FiPy| is an object oriented, partial differential equation (PDE) solver,
written in Python_ |citePython|, based on a standard finite volume (FV) approach.  The
framework has been developed in the `Metallurgy Division`_ and Center for
Theoretical and Computational Materials Science (CTCMS_), in the Materials
Science and Engineering Laboratory (MSEL_) at the National Institute of
Standards and Technology (NIST_).

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
coupled elliptic, hyperbolic and parabolic PDEs.  Current models include
phase field |citePhaseField| treatments of polycrystalline, dendritic, and electrochemical
phase transformations as well as a level set treatment of the
electrodeposition process |citeCEAC|.

The latest information about |FiPy| can be found at
http://www.ctcms.nist.gov/fipy/.

--------
Download
--------

|FiPy| is available for download via CVS_ or as a `compressed
archive`_.  |FiPy| can be redistributed and/or modified freely, provided
that any derivative works bear some notice that they are derived from it,
and any modified versions bear some notice that they have been modified.

------------
Installation
------------

Please refer to |INSTALLATION.txt|.

-------
Support
-------

|FiPy| is being actively developed and supported.  Please use the `tracking
system`_ for bugs, support requests, feature requests and patch
submissions.  A `discussion forum`_ is also available.  We are also seeking
collaborative efforts on this project.

------------------------
Conventions and Notation
------------------------

|FiPy| is driven by Python_ script files than you can view or modify in any
text editor.  |FiPy| sessions are invoked from a command-line shell, such
as ``tcsh`` or ``bash``.

Throughout, text to by typed at the keyboard will appear ``like this``.
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

   Python_ is sensitive to indentation and care should be take to enter
   text exactly as it appears in the examples.

When references are made to file system paths, it is assumed that the
current working directory is the |FiPy| distribution directory, refered to
as the "base directory", such that::

    examples/diffusion/steadyState/mesh1D/input.py

will correspond to, *e.g.*::

    /some/where/FiPy-0.1/examples/diffusion/steadyState/mesh1D/input.py

Paths will always be rendered using POSIX conventions.  Any references of
the form::

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
.. _Python:               http://www.python.org/
.. _CVS:                  http://cvs.sourceforge.net/viewcvs.py/fipy/
.. _compressed archive:   http://sourceforge.net/project/showfiles.php?group_id=118428
.. _tracking system:      http://sourceforge.net/tracker/?group_id=118428&atid=681142
.. _discussion forum:     http://sourceforge.net/forum/?group_id=118428

.. include:: utils/include.txt

.. |FiPy| replace:: |htmlFiPy| |latexFiPy|
.. |INSTALLATION.txt| replace:: |latexINSTALLATION.txt| |htmlINSTALLATION.txt|


.. |citePython| raw:: latex

   \cite{Python}

.. |citePhaseField| raw:: latex

   \cite{BoettingerReview:2002,McFaddenReview:2002}

.. |citeCEAC| raw:: latex

   \cite{NIST:damascene:2001}
