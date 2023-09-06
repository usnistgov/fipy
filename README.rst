========
Overview
========



| |Tests|_ |Documentation|_ |nix|_
| |GitHub|_ |PyPI|_  |CondaForge|_ |Binder|_
| |gitter|_ |OpenHub|_

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

   See the latest updates in the :ref:`CHANGELOG`.

---------------------------------
Even if you don't read manuals...
---------------------------------

...please read :ref:`INSTALLATION`, :ref:`USAGE` and :ref:`FAQ`, as well
as :mod:`examples.diffusion.mesh1D`.

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

We offer several modes to communicate with the :term:`FiPy` developers and
with other users.

.. toctree::

   CONTACT

We welcome collaborative efforts on this project.

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
current working directory is the :term:`FiPy` distribution directory, referred to
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

.. |GitHub|        image:: https://img.shields.io/github/contributors/usnistgov/fipy.svg
.. _GitHub:        https://github.com/usnistgov/fipy
.. |gitter|        image:: https://badges.gitter.im/usnistgov/fipy.svg
.. _gitter:        https://gitter.im/usnistgov/fipy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=body_badge
.. |Tests|         image:: https://dev.azure.com/guyer/FiPy/_apis/build/status/usnistgov.fipy?branchName=master
.. _Tests:         https://dev.azure.com/guyer/FiPy/_build?definitionId=2
.. |Documentation| image:: https://github.com/usnistgov/fipy/actions/workflows/NISTtheDocs2Death.yml/badge.svg
.. _Documentation: https://github.com/usnistgov/fipy/actions/workflows/NISTtheDocs2Death.yml
.. |nix|           image:: https://github.com/usnistgov/fipy/actions/workflows/nix.yml/badge.svg
.. _nix:           https://github.com/usnistgov/fipy/actions/workflows/nix.yml
.. |OpenHub|       image:: https://www.openhub.net/p/fipy/widgets/project_thin_badge.gif
.. _OpenHub:       https://www.openhub.net/p/fipy
.. |PyPI|          image:: https://img.shields.io/pypi/v/fipy.svg
.. _PyPI:          https://pypi.python.org/pypi/FiPy
.. |CondaForge|    image:: https://img.shields.io/conda/pn/conda-forge/fipy?label=conda-forge
.. _CondaForge:    https://anaconda.org/conda-forge/fipy

.. |Binder|        image:: https://mybinder.org/badge.svg
.. _Binder:        https://mybinder.org/v2/gh/usnistgov/fipy/master?filepath=examples%2Findex.ipynb
