========
Overview
========

|FiPy| is an object oriented partial differential equation (PDE) solver,
written in Python_, based on a standard finite volume (FV) approach.  The
initial development of the framework has been undertaken in the Center for
Theoretical and Computational Materials Science (CTCMS_), in the Materials
Science and Engineering Laboratory (MSEL_) at the National Institute of
Standards and Technology (NIST_).

The solution of coupled sets of PDEs is ubiquitous in the numerical
simulation of science problems.  Numerous PDE solvers exist using a
variety of languages and numerical approaches. Many are proprietary,
expensive and difficult to customize.  As a result, scientists spend
considerable resources repeatedly developing limited tools for
specific problems.  Our approach, combining the FV method and Python_,
provides a tool that is extensible, powerful and freely available. A
significant advantage to Python_ is the existing suite of tools for
array calculations, sparse matrices and data representation.

The |FiPy| framework includes terms for transient diffusion, convection and
standard sources, enabling the solution of arbitrary combinations of
coupled elliptic, hyperbolic and parabolic PDEs.  Current models include
phase field treatments of polycrystalline, dendritic, and electrochemical
phase transformations as well as a level set treatment of the
electrodeposition process.

--------
Download
--------

|FiPy| is freely available for download from sourceforge_ via CVS_ or as
a `compressed archive`_. |FiPy| can be redistributed and/or modified freely provided
that any derivative works bear some notice that they are derived from
it, and any modified versions bear some notice that they have been
modified.

------------
Installation
------------

Please refer to |INSTALLATION.txt|.

-------
Support
-------

The developers are committed to supporting the project. The project is
hosted via sourceforge_. Please use the sourceforge_ `tracking system`_
for bugs, support requests, feature requests and patch submissions. A
`discussion forum`_ is also available. We are also seeking collaborative
efforts on this project. Please get involved.

.. _MSEL:                 http://www.msel.nist.gov/
.. _CTCMS:                http://www.ctcms.nist.gov/
.. _NIST:                 http://www.nist.gov/
.. _Python:               http://www.python.org/
.. _sourceforge:          http://sourceforge.net/projects/fipy/
.. _CVS:                  http://cvs.sourceforge.net/viewcvs.py/fipy/
.. _compressed archive:   http://sourceforge.net/project/showfiles.php?group_id=118428
.. _tracking system:      http://sourceforge.net/tracker/?group_id=118428&atid=681142
.. _discussion forum:     http://sourceforge.net/forum/?group_id=118428

.. include:: documentation/tools/include.txt

.. |FiPy| replace:: |latexFiPy| |htmlFiPy|
.. |INSTALLATION.txt| replace:: |latexINSTALLATION.txt| |htmlINSTALLATION.txt|

