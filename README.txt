--------
Overview
--------

|FiPy| is an object oriented partial differential equation (PDE)
solver written in Python_ based on a standard finite volume (FV)
approach.

The solution of coupled sets of PDEs is ubiquitous in the numerical
simulation of science problems.  Numerous PDE solvers exist using a
variety of languages and numerical approaches. Many are proprietary,
expensive and difficult to customize.  As a result, scientists spend
considerable resources repeatedly developing limited tools for
specific problems.  Our approach, combining the FV method and Python_,
provides a tool that is extensible, powerful and freely available. A
significant advantage to Python_ is the existing suite of tools for
array calculations, sparse matrices and data representation.

The |FiPy| framework includes terms for transient diffusion,
convection and standard sources, enabling the solution of arbitrary
combinations of coupled elliptic, hyperbolic and parabolic
PDEs. Current models include phase field treatments of
electrochemical, polycrystalline and dendritic phase transformations
as well as a level set treatment of the electrodeposition process.

----
Aims
----

The goal of the |FiPy| framework is to develop a highly customizable
open source code available to the scientific community. |FiPy| allows
users to select and customize modules from within the framework. Those
modules can then be combined with an extremely powerful high level
scripting language allowing the user full control over the solution
algorithm. This 'top down' method of computer programming allows
significant code reuse while accounting for the unique vagaries of
each problem.

The initial development of the framework has been undertaken in the
Center for Theoretical and Computational Materials Science (CTCMS_),
in the Materials Science and Engineering Laboratory (MSEL_) at the
National Institute of Standards and Technology (NIST_). The framework
is known as |FiPy| and is developed at the top level in the Python_
programming language. Python_ enables the integration of low-level
high performance languages such as C and FORTRAN allowing for
efficient optimization. Python_ also has a large repository of free
numerical and scientific software that can be integrated into the
|FiPy| framework. At NIST_, there is a significant ongoing effort to
maximize the efficiency of |FiPy| with the use of third party software
and the judicious use of the C programming language embedded in
standard core modules.

The end goal of |FiPy| is to provide a realistic open source choice
for solving coupled sets of equations. Once an efficient customizable
framework is developed, the end goal will only be achieved with proper
support for the community of users. Many previous scientific computing
projects, although valuable tools, are not adopted by the end users
due to inadequate distribution tools and lack of support. Adequate
support consists of clear documentation, efficient distribution
methods, a code repository, regular code updates (bug fixes) and a
test infrastructure. We aim to address all of these issues.

The aims mentioned above can be summarized as follows:

* Develop core test piece models for various material science problems.
* Develop a highly efficient and customizable code.
* Support the distribution of code to the end users.

.. |FiPy| raw:: latex

	  \FiPy{}

.. _Python:   http://www.python.org/
.. _MSEL:     http://www.msel.nist.gov/
.. _CTCMS:    http://www.ctcms.nist.gov/
.. _NIST:     http://www.nist.gov/
