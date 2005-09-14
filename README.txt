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

-------------------------
What's new in version 1.0
-------------------------

Numerous changes have been made since |FiPy| 0.1 was released, but the most 
signficant ones are:

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

- A new viewer based on `matplotlib`_ has been added.

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

The syntax of |FiPy| 1.0 scripts is incompatible with earlier releases.  A 
tutorial for updating your existing scripts can be found in 
|examples/update0_1to1_0.py|.

.. _matplotlib:  http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://matplotlib.sourceforge.net

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

|FiPy| is being actively developed and supported.  Please use the `tracking
system`_ for bugs, support requests, feature requests and patch
submissions.  A `mailing list`_ is also available.  We are also seeking
collaborative efforts on this project.

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
.. _Python:               http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://www.python.org/
.. _CVS:                  http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://cvs.sourceforge.net/viewcvs.py/fipy/
.. _compressed archive:   http://www.ctcms.nist.gov/fipy/download/FiPy-0.1.tar.gz
.. _tracking system:      http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://sourceforge.net/tracker/?group_id=118428
.. _mailing list:         http://www.ctcms.nist.gov/fipy/mail.html

.. include:: utils/include.txt

.. |FiPy| replace:: |htmlFiPy| |latexFiPy|
.. |INSTALLATION-txt| replace:: |htmlINSTALL| |latexINSTALL|
.. |the FAQ| replace:: |htmlFAQ| |latexFAQ|

.. |latexUpdate0_1to1_0.py| raw:: latex

   Chapter~\ref{chap:Update0.1to1.0}

.. |htmlUpdate0_1to1_0.py| raw:: html

   <a href="http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://cvs.sourceforge.net/viewcvs.py/fipy/fipy/examples/update0_1to1_0.py?view=markup">examples/update0_1to1_0.py</a>

.. |examples/update0_1to1_0.py| replace:: |latexUpdate0_1to1_0.py| |htmlUpdate0_1to1_0.py|

.. |citePython| raw:: latex

   \cite{Python}

.. |citePhaseField| raw:: latex

   \cite{BoettingerReview:2002,McFaddenReview:2002}

.. |citeCEAC| raw:: latex

   \cite{NIST:damascene:2001}
