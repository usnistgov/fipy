.. image:: https://raw.githubusercontent.com/usnistgov/fipy/master/docs/source/_static/logo.png
  :target: https://pages.nist.gov/fipy/en/stable
  :width: 110
  :height: 110
  :align: left

| |Tests|_ |Documentation|_ |nix|_
| |GitHub|_ |PyPI|_  |CondaForge|_ |Binder|_
| |gitter|_ |OpenHub|_

- **Documentation:** https://pages.nist.gov/fipy/en/stable
- **Development version of the documentation:** https://pages.nist.gov/fipy/en/latest
- **Source code:** https://github.com/usnistgov/fipy
- **Installation:** https://pages.nist.gov/fipy/en/latest/INSTALLATION.html
- **Usage:** https://pages.nist.gov/fipy/en/latest/USAGE.html
- **Examples:** https://pages.nist.gov/fipy/en/stable/EXAMPLES.html
- **Support:** https://pages.nist.gov/fipy/en/stable/README.html#support
- **FAQ:** https://pages.nist.gov/fipy/en/latest/FAQ.html
- **Releases:** https://pages.nist.gov/fipy/en/stable/CHANGELOG.html

FiPy is an object oriented, partial differential equation (PDE)
solver, written in Python, based on a standard finite volume
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
specific problems.  Our approach, combining the FV method and Python,
provides a tool that is extensible, powerful and freely available. A
significant advantage to Python is the existing suite of tools for
array calculations, sparse matrices and data rendering.

The FiPy framework includes terms for transient diffusion, convection and
standard sources, enabling the solution of arbitrary combinations of
coupled elliptic, hyperbolic and parabolic PDEs.  Currently implemented
models include phase field treatments of polycrystalline, dendritic, and
electrochemical phase transformations, as well as drug eluting stents,
reactive wetting.  photovoltaics and a level set treatment of the
electrodeposition process.

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
.. |Documentation| image:: https://github.com/usnistgov/fipy/actions/workflows/Docs4NIST.yml/badge.svg
.. _Documentation: https://github.com/usnistgov/fipy/actions/workflows/Docs4NIST.yml
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

