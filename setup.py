""" FiPy is an object oriented, partial differential equation (PDE) solver

FiPy is based on a standard finite volume (FV) approach.  The framework has
been developed in the Materials Science and Engineering Division (MSED) and
Center for Theoretical and Computational Materials Science (CTCMS), in the
Material Measurement Laboratory (MML) at the National Institute of
Standards and Technology (NIST).
"""

from setuptools import setup

import versioneer

VERSION = versioneer.get_version()

DIST = setup(
    version=VERSION,
    cmdclass=dict(
        **versioneer.get_cmdclass()
    ),
)
