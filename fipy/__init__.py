"""
:term:`FiPy` is an object oriented, partial differential equation (PDE) solver,
written in :term:`Python`, based on a standard finite volume (FV) approach. The
framework has been developed in the Materials Science and Engineering Division
(MSED_) and Center for Theoretical and Computational Materials Science (CTCMS_),
in the Material Measurement Laboratory (MML_) at the National Institute of
Standards and Technology (NIST_).

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

.. _MML:                  http://www.nist.gov/mml/
.. _CTCMS:                http://www.ctcms.nist.gov/
.. _MSED:                 http://www.nist.gov/mml/msed/
.. _NIST:                 http://www.nist.gov/
"""
from __future__ import unicode_literals
from builtins import input
__docformat__ = 'restructuredtext'

import sys

from fipy.boundaryConditions import *
from fipy.meshes import *
from fipy.solvers import *
from fipy.steppers import *
from fipy.terms import *
from fipy.tools import *
from fipy.variables import *
from fipy.viewers import *

__all__ = []
__all__.extend(boundaryConditions.__all__)
__all__.extend(meshes.__all__)
__all__.extend(solvers.__all__)
__all__.extend(steppers.__all__)
__all__.extend(terms.__all__)
__all__.extend(tools.__all__)
__all__.extend(variables.__all__)
__all__.extend(viewers.__all__)

# fipy needs to export raw_input whether or not parallel

input_original = input

if parallelComm.Nproc > 1:
    def mpi_input(prompt=""):
        parallelComm.Barrier()
        sys.stdout.flush()
        if parallelComm.procID == 0:
            sys.stdout.write(prompt)
            sys.stdout.flush()
            return sys.stdin.readline()
        else:
            return ""
    input = mpi_input

__all__.extend(['input', 'input_original'])

from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

_saved_stdout = sys.stdout

def _serial_doctest_raw_input(prompt):
    """Replacement for `raw_input()` that works in doctests
    """
    _saved_stdout.write("\n")
    _saved_stdout.write(prompt)
    _saved_stdout.flush()
    return sys.stdin.readline()

def doctest_raw_input(prompt):
    """Replacement for `raw_input()` that works in doctests

    This routine attempts to be savvy about running in parallel.
    """
    try:
        from fipy.tools import parallelComm
        parallelComm.Barrier()
        _saved_stdout.flush()
        if parallelComm.procID == 0:
            txt = _serial_doctest_raw_input(prompt)
        else:
            txt = ""
        parallelComm.Barrier()
    except ImportError:
        txt = _serial_doctest_raw_input(prompt)
#     return txt

def test(*args):
    r"""
    Test `Fipy`. Equivalent to::

        $ python setup.py test --modules

    Use

    >>> import fipy
    >>> fipy.test('--help')

    for a full list of options. Options can be passed in the same way
    as they are appended at the command line. For example, to test
    `FiPy` with `Trilinos` and inlining switched on, use

    >>> fipy.test('--trilinos', '--inline')

    At the command line this would be::

        $ python setup.py test --modules --trilinos --inline

    """

    from setuptools import setup
    from fipy.tests.test import test
    import tempfile

    tmpDir = tempfile.mkdtemp()

    try:
        setup(name='FiPy',
              script_args = ['egg_info', '--egg-base=' + tmpDir,
                             'test', '--modules'] + list(args),
              cmdclass={'test': test})
    except SystemExit as exitErr:
        import shutil
        shutil.rmtree(tmpDir)
        raise exitErr

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
