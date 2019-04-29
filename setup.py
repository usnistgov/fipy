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

from _setup.build_docs import build_docs
from _setup.changelog import changelog
from _setup.copy_script import copy_script
from _setup.upload_products import upload_products
from _setup.release import release

# bootstrap setuptools for users that don't already have it
import ez_setup

ez_setup.use_setuptools()

try:
    FILE = open("README.rst", "r")
    LONG_DESCRIPTION = "\n" + FILE.read() + "\n"
    FILE.close()
except IOError as _:
    LONG_DESCRIPTION = ""

try:
    FILE = open("LICENSE.rst", "r")
    LICENSE = "\n" + "".join([" " * 8 + l for l in FILE])
    FILE.close()
except IOError as _:
    LICENSE = ""

DIST = setup(
    name="FiPy",
    install_requires=["numpy", "scipy", "matplotlib"],
    version=versioneer.get_version(),
    author="Jonathan Guyer, Daniel Wheeler, & Jim Warren",
    author_email="fipy@nist.gov",
    url="http://www.ctcms.nist.gov/fipy/",
    license=LICENSE,
    description="A finite volume PDE solver in Python",
    long_description=LONG_DESCRIPTION,
    cmdclass=dict(
        {
            "build_docs": build_docs,
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
