import os

from distutils.core import Command
import versioneer

from fipy.tools.performance.efficiency_test import Efficiency_test
from fipy.tools.copy_script import Copy_script
from fipy.tools.changelog import changelog
from fipy.tests.testClass import _TestClass

# bootstrap setuptools for users that don't already have it
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

# from fipy.tests.testRunner import test, unittest

from setuptools.command.test import test as _test

            
test = _TestClass(_test)

try:
    # we only need "unittest" if bitten is installed 
    # (and we're running as a bitten.slave)
    from bitten.util.testrunner import unittest as _unittest
    unittest = _TestClass(_unittest)
except ImportError, e:
    unittest = test



class build_docs(Command):

    description = "build the FiPy documentation"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [('pdf', None, "compile the PDF variant of the documentation"),
                    ('html', None, "compile the HTML variant of the documentation"),
                    ('cathartic', None, "rewrite all the files (default is to only rewrite changed files)"),
                   ]

    def initialize_options (self):
        self.pdf = 0
        self.html = 0
        self.cathartic = 0

    def finalize_options (self):
        pass

    def run (self):
        import sphinx.cmd.build
        import sphinx.ext.apidoc
        
        sphinx_args = ['-P', '-n', '-c', 'documentation/', '.']
        apidoc_args = []
        
        if self.cathartic:
            sphinx_args = ['-a', '-E'] + sphinx_args
            apidoc_args = ['--force'] + apidoc_args
            
        sphinx.ext.apidoc.main(['--output-dir=fipy/generated', '--suffix=rst']
                    + apidoc_args + ['fipy'])
        sphinx.ext.apidoc.main(['--output-dir=documentation/tutorial/package/generated', '--suffix=rst']
                    + apidoc_args + ['documentation/tutorial/package'])

        if self.html:
            sphinx.cmd.build.main(['-b', 'redirecting_html'] + sphinx_args + ['documentation/_build/html/'])

        if self.pdf:
            try:
                sphinx.cmd.build.main(['-b', 'latex'] + sphinx_args + ['documentation/_build/latex/'])
            except SystemExit:
                pass
            
            outdir = os.path.join('documentation', '_build', 'latex')
            
            from docutils.core import publish_file

            for xtra in ("LICENSE", "DISCLAIMER"):
                publish_file(source_path="%s.rst" % xtra,
                             destination_path=os.path.join(outdir, "%s.tex" % xtra),
                             reader_name='standalone',
                             parser_name='restructuredtext',
                             writer_name='latex',
                             settings_overrides= {
                                 'template': 'documentation/_templates/empty.tex'
                             })

            savedir = os.getcwd()
            
            os.chdir(outdir)
                
            os.system("pdflatex fipy")
            os.system("pdflatex fipy")
            os.system("pdflatex fipy")
            os.system("makeindex -s python.ist fipy")
            os.system("makeindex -s python.ist modfipy")
            os.system("pdflatex fipy")
            os.system("pdflatex fipy")
                
            os.chdir(savedir)
            
class upload_products(Command):
    description = "upload FiPy compressed archives to website(s)"
    
    user_options = [('pdf', None, "upload the PDF variant of the documentation"),
                    ('html', None, "upload the HTML variant of the documentation"),
                    ('tarball', None, "upload the .tar.gz source distribution"),
                    ('winzip', None, "upload the .win32.zip distribution"),
                   ]

    def initialize_options (self):
        self.pdf = 0
        self.html = 0
        self.tarball = 0
        self.winzip = 0

    def finalize_options (self):
        pass

    def run(self):
        if self.pdf:
            print "setting permissions of manual..."
            os.system('chmod -R g+w documentation/_build/latex/fipy.pdf')
            
            print "linking manual to `dist/`..."
            os.system('mkdir dist/')
            os.system('ln -f documentation/_build/latex/fipy.pdf dist/fipy-%s.pdf'%self.distribution.metadata.get_version())
            
        if self.html:
            print "setting group and ownership of web pages..."
            os.system('chmod -R g+w documentation/_build/html/')
            
            print "uploading web pages..."
            # The -t flag (implicit in -a) is suddenly causing problems
            # os.system('rsync -aLC -e ssh %s %s'%('documentation/www/', os.environ['FIPY_WWWHOST']))
            os.system('rsync -rlpgoDLC -e ssh %s %s' % ('documentation/_build/html/', os.environ['FIPY_WWWHOST']))

        if self.tarball:
            file = 'dist/FiPy-%s.tar.gz' % self.distribution.metadata.get_version()
            print "setting permissions for %s ..." % file
            os.system('chmod -R g+w %s' % file)

            print "uploading tarball..."
            os.system('rsync -pgoDLC -e ssh %s %s/download/' % (file, os.environ['FIPY_WWWHOST']))

        if self.winzip:
            file = 'dist/FiPy-%s.win32.zip' % self.distribution.metadata.get_version()
            print "setting permissions for %s ..." % file
            os.system('chmod -R g+w %s' % file)
            
            print "uploading winzip..."
            os.system('rsync -pgoDLC -e ssh %s %s/download/' % (file, os.environ['FIPY_WWWHOST']))

        if self.pdf or self.tarball or self.winzip:
            print "activating web pages..."
            os.system(os.environ['FIPY_WWWACTIVATE'])

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
            'build_docs':build_docs,
            'upload_products':upload_products,
            'test':test,
            'unittest':unittest,
            'copy_script': Copy_script,
            'efficiency_test': Efficiency_test,
            'changelog': changelog,
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
