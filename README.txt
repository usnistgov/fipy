--------
Overview
--------

FiPy is an object oriented partial differential equation (PDE) solver
written in Python_ based on a standard finite volume (FV) approach.

The solution of coupled sets of PDEs is ubiquitous in the numerical
simulation of science problems.  Numerous PDE solvers exist using a
variety of languages and numerical approaches. Many are proprietary,
expensive and difficult to customize.  As a result, scientists spend
considerable resources repeatedly developing limited tools for
specific problems.  Our approach, combining the FV method and Python_,
provides a tool that is extensible, powerful and freely available. A
significant advantage to Python_ is the existing suite of tools for
array calculations, sparse matrices and data representation.

The FiPy framework includes terms for transient diffusion, convection
and standard sources, enabling the solution of arbitrary combinations
of coupled elliptic, hyperbolic and parabolic PDEs. Current models
include phase field treatments of electrochemical, polycrystalline and
dendritic phase transformations as well as a level set treatment of
the electrodeposition process.

--------
Download
--------

FiPy is freely available for download from sourceforge_ via CVS_ or as
a tarball.

------------
Installation
------------

Please refer to the installation guide_ or the manual_.

------
Manual
------

The manual_ can be downloaded and is also part of the distribution
(documentation/manual/fipy.pdf).

A fresh copy of the manual_ can be built with the use of epydoc_ by
issuing the following commands in the base directory::

     python setup.py build_docs --latex --manual
     python setup.py build_docs --manual

-------
Support
-------

The developers are committed to supporting the project. The project is
hosted via sourceforge_. Please use the sourceforge_ tracking system_
for bugs, support requests, feature requests and patch submissions. A
discussion forum_ is also available.

.. |FiPy| raw:: latex

	  \FiPy{}

.. _Python:      http://www.python.org/
.. _MSEL:        http://www.msel.nist.gov/
.. _CTCMS:       http://www.ctcms.nist.gov/
.. _NIST:        http://www.nist.gov/
.. _webpage:     index.html
.. _CVS:         http://cvs.sourceforge.net/viewcvs.py/fipy/
.. _sourceforge: http://sourceforge.net/projects/fipy/
.. _manual:      ../manual/fipy.pdf
.. _guide:       installation.html
.. _bug tracker: http://sourceforge.net/tracker/?group_id=118428&atid=681141
.. _system:      http://sourceforge.net/tracker/?group_id=118428&atid=681142
.. _forum:       http://sourceforge.net/forum/?group_id=118428
.. _epydoc:      http://epydoc.sourceforge.net/