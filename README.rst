.. |.changelog| replace:: Change Log
.. _.changelog: https://pages.nist.gov/fipy/en/latest/CHANGELOG.html#changelog
.. |.examples.diffusion.mesh1D| replace:: ``examples.diffusion.mesh1D``
.. _.examples.diffusion.mesh1D: https://github.com/usnistgov/fipy/blob/f0ee6b25fb5f1e7063be3a5195d547c83cfd0ddc/examples/diffusion/mesh1D.py
.. |.faq| replace:: Frequently Asked Questions
.. _.faq: https://pages.nist.gov/fipy/en/latest/FAQ.html#faq
.. |.FiPy| replace:: FiPy
.. _.FiPy: https://pages.nist.gov/fipy/en/latest/glossary.html#term-FiPy
.. |.installation| replace:: Installation
.. _.installation: https://pages.nist.gov/fipy/en/latest/INSTALLATION.html#installation
.. |.Python| replace:: Python
.. _.Python: https://pages.nist.gov/fipy/en/latest/glossary.html#term-Python
.. |.usage| replace:: Using FiPy
.. _.usage: https://pages.nist.gov/fipy/en/latest/USAGE.html#usage


========
Overview
========



| |Tests|_ |Documentation|_ |nix|_
| |GitHub|_ |PyPI|_  |CondaForge|_ |Binder|_
| |gitter|_ |OpenHub|_

|.FiPy|_ is an object oriented, partial differential equation (PDE)
solver, written in |.Python|_, based on a standard finite volume
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
specific problems.  Our approach, combining the FV method and |.Python|_,
provides a tool that is extensible, powerful and freely available. A
significant advantage to |.Python|_ is the existing suite of tools for
array calculations, sparse matrices and data rendering.

The |.FiPy|_ framework includes terms for transient diffusion,
convection and standard sources, enabling the solution of arbitrary
combinations of coupled elliptic, hyperbolic and parabolic PDEs. Currently
implemented models include phase field :cite:`BoettingerReview:2002`
:cite:`ChenReview:2002` :cite:`McFaddenReview:2002` treatments of polycrystalline,
dendritic, and electrochemical phase transformations, as well as drug
eluting stents :cite:`Saylor:2011p2794`, reactive wetting :cite:`PhysRevE.82.051601`,
photovoltaics :cite:`Hangarter:2011p2795` and a level set treatment of the
electrodeposition process :cite:`NIST:damascene:2001`.

The latest information about |.FiPy|_ can be found at
http://www.ctcms.nist.gov/fipy/.

See the latest updates in the |.changelog|_.

---------------------------------
Even if you don't read manuals...
---------------------------------

...please read |.installation|_, |.usage|_ and |.faq|_, as well
as |.examples.diffusion.mesh1D|_.

-------------------------
Download and Installation
-------------------------

Please refer to |.installation|_ for details on download and
installation. |.FiPy|_ can be redistributed and/or modified
freely, provided that any derivative works bear some notice that they
are derived from it, and any modified versions bear some notice that
they have been modified.

-------
Support
-------

We offer several modes to communicate with the |.FiPy|_ developers and
with other users.

* `Contact <https://pages.nist.gov/fipy/en/latest/CONTACT.html>`_

  * `GitHub Discussions <https://pages.nist.gov/fipy/en/latest/CONTACT.html#github-discussions>`_
  * `GitHub Issues <https://pages.nist.gov/fipy/en/latest/CONTACT.html#github-issues>`_
  * `StackOverflow <https://pages.nist.gov/fipy/en/latest/CONTACT.html#stackoverflow>`_
  * `Mailing List <https://pages.nist.gov/fipy/en/latest/CONTACT.html#mailing-list>`_

    * `List Archive <https://pages.nist.gov/fipy/en/latest/CONTACT.html#list-archive>`_



      |



We welcome collaborative efforts on this project.

------------------------
Conventions and Notation
------------------------

|.FiPy|_ is driven by |.Python|_ script files than you can view or modify in any
text editor.  |.FiPy|_ sessions are invoked from a command-line shell, such
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

are intended to indicate an interactive session in the |.Python|_ interpreter.
We will refer to these as "interactive sessions" or as "doctest blocks".
The text "``>>>``" at the beginning of a line denotes the *primary prompt*,
calling for input of a |.Python|_ command.  The text "``...``" denotes the
*secondary prompt*, which calls for input that continues from the line
above, when required by |.Python|_ syntax.  All remaining lines, which begin
at the left margin, denote output from the |.Python|_ interpreter.  In all
cases, the prompt is supplied by the |.Python|_ interpreter and should not be
typed by you.


.. list-table::
   :header-rows: 1
   
   * - 🚩 Warning
   * - |.Python|_ is sensitive to indentation and care should be taken to enter
       text exactly as it appears in the examples.


When references are made to file system paths, it is assumed that the
current working directory is the |.FiPy|_ distribution directory, referred to
as the "base directory", such that::

    examples/diffusion/steadyState/mesh1D.py

will correspond to, *e.g.*::

    /some/where/FiPy-X.Y/examples/diffusion/steadyState/mesh1D.py

Paths will always be rendered using POSIX conventions (path elements
separated by "``/``").  Any references of the form::

    examples.diffusion.steadyState.mesh1D

are in the |.Python|_ module notation and correspond to the equivalent POSIX
path given above.

We may at times use a


.. list-table::
   :header-rows: 1
   
   * - 📝 Note
   * - to indicate something that may be of interest


or a


.. list-table::
   :header-rows: 1
   
   * - 🚩 Warning
   * - to indicate something that could cause serious problems.


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

