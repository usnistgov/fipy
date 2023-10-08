""" FiPy is an object oriented, partial differential equation (PDE) solver

FiPy is based on a standard finite volume (FV) approach.  The framework has
been developed in the Materials Science and Engineering Division (MSED) and
Center for Theoretical and Computational Materials Science (CTCMS), in the
Material Measurement Laboratory (MML) at the National Institute of
Standards and Technology (NIST).
"""
from __future__ import unicode_literals

from setuptools import setup, find_packages

import versioneer

from _setup.changelog import changelog
from _setup.copy_script import copy_script
from _setup.upload_products import upload_products
from _setup.release import release

LONG_DESCRIPTION = """
FiPy is an object oriented, partial differential equation (PDE) solver,
written in Python, based on a standard finite volume (FV) approach.  This
combination provides a tool that is extensible, powerful and freely
available.  A significant advantage to Python is the existing suite of
tools for array calculations, sparse matrices and data rendering.

The FiPy framework includes terms for transient diffusion, convection and
standard sources, enabling the solution of arbitrary combinations of
coupled elliptic, hyperbolic and parabolic PDEs.  Currently implemented
models include phase field treatments of polycrystalline, dendritic, and
electrochemical phase transformations, as well as drug eluting stents,
reactive wetting, photovoltaics and a level set treatment of the
electrodeposition process.
"""

VERSION = versioneer.get_version()

DIST = setup(
    name="FiPy",
    install_requires=["numpy", "scipy", "matplotlib", "future"],
    version=VERSION,
    download_url="https://github.com/usnistgov/fipy/archive/{}.zip".format(VERSION),
    author="Jonathan Guyer, Daniel Wheeler, & Jim Warren",
    author_email="fipy@list.nist.gov",
    url="http://www.ctcms.nist.gov/fipy/",
    license="NIST Public Domain",
    description="A finite volume PDE solver in Python",
    long_description=LONG_DESCRIPTION,
    cmdclass=dict(
        {
            "upload_products": upload_products,
            "copy_script": copy_script,
            "changelog": changelog,
            "release": release,
        },
        **versioneer.get_cmdclass()
    ),
    test_suite="fipy.testFiPy._suite",
    packages=find_packages(exclude=["examples", "examples.*", "utils", "utils.*"]),
    entry_points="""
                 [fipy.viewers]
                 matplotlib = fipy.viewers.matplotlibViewer:MatplotlibViewer
                 mayavi = fipy.viewers.mayaviViewer:MayaviClient
                 [distutils.commands]
                 test = fipy.tests.test:test
                 unittest = fipy.tests.test:test
             """,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Environment :: X11 Applications",
        "Intended Audience :: Science/Research",
        "License :: Public Domain",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
