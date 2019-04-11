import os

from distutils.core import Command
import versioneer

from _setup.build_docs import build_docs
from _setup.changelog import changelog
from _setup.copy_script import copy_script
from _setup.testClass import _TestClass
from _setup.upload_products import upload_products

# bootstrap setuptools for users that don't already have it
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

from setuptools.command.test import test as _test

            
test = _TestClass(_test)


try:            
    f = open('README.rst', 'r')
    long_description = '\n' + f.read() + '\n'
    f.close()
except IOError, e:
    long_description = ''
        
try:
    f = open('LICENSE.rst', 'r') 
    license = '\n' + ''.join([' '*8 + l for l in f])
    f.close()
except IOError, e:
    license = ''    
# The following doesn't work reliably, because it requires fipy
# to already be installed (or at least egged), which is kind of 
# obnoxious. We use cmdclass instead.
# 
#         entry_points = {
#             'distutils.commands': [
#                 'test = fipy.tests.testRunner:test',
#                 'unittest = fipy.tests.testRunner:unittest', 
#             ],
#         },

##Hacked from numpy
def getVersion():
    import subprocess
    import warnings

    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'

        out = subprocess.check_output(cmd, env=env)
        # ticket:475 - fix for bytecode received in Py3k
        # http://jeetworks.org/node/67
        out = out.decode("utf-8")
        return out

    version = 'unknown'

    if os.path.exists('.git'):
        try:
            out = _minimal_ext_cmd(['git', 'describe', '--tags', '--match', 'version-*'])
            # convert git long-form version string, e.g., "version-3_1_1-127-g413ed61",
            # into PEP 440 version, e.g., "3.1.1.dev127+g413ed61"
            version = out.strip().split("-")
            suffix = version[2:]
            version = ".".join(version[1].split("_"))
            if suffix:
                dev, sha = suffix
                version = "%s.dev%s+%s" % (version, dev, sha)
        except OSError:
            import warnings
            warnings.warn("Could not run ``git describe``")
        except subprocess.CalledProcessError:
            warnings.warn("Could not read git tags")
            version = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
            version = version.strip()

    if (version == 'unknown') and os.path.exists('FiPy.egg-info'):
        from fipy import _getVersion
        version = _getVersion()

    return version

dist = setup(	name = "FiPy",
        version = versioneer.get_version(),
        author = "Jonathan Guyer, Daniel Wheeler, & Jim Warren",
        author_email = "fipy@nist.gov",
        url = "http://www.ctcms.nist.gov/fipy/",
        license = license,
        description = "A finite volume PDE solver in Python",
        long_description = long_description,
        cmdclass = dict({
            'build_docs': build_docs,
            'upload_products': upload_products,
            'test': test,
            'unittest': test,
            'copy_script': copy_script,
            'changelog': changelog
            },
            **versioneer.get_cmdclass()
        ),
        test_suite="fipy.testFiPy._suite",
        packages = find_packages(exclude=["examples", "examples.*", "utils", "utils.*"]),
        entry_points="""
            [fipy.viewers]
            matplotlib = fipy.viewers.matplotlibViewer:MatplotlibViewer
            mayavi = fipy.viewers.mayaviViewer:MayaviClient
        """,
        classifiers = [
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Environment :: X11 Applications',
            'Intended Audience :: Science/Research',
            'License :: Public Domain',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Scientific/Engineering :: Visualization',
            'Topic :: Software Development :: Libraries :: Python Modules'
        ],
)

if 'install' in dist.commands:
    req = []
    
    for pkg in ['numpy', 'pysparse']:
        try:
            __import__(pkg)
        except ImportError, exc:
            req.append(pkg)
            
    if len(req) > 0:
        print "!!!!!!"
        print "The required module(s) " + str(req) + " cannot be loaded."
        print "FiPy will not work properly until these modules are installed."

    opt = []
    
    for pkg in ['scipy', 'matplotlib', 'mayavi']:
        try:
            __import__(pkg)
        except ImportError, exc:
            opt.append(pkg)
        
    if len(opt) > 0:
        print "------"
        print "The optional module(s) " + str(opt) + " cannot be loaded."
        print "FiPy will have improved capabilities if these modules are installed."
