#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "setup.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  setup.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##


import glob
import os
import sys
import string

from distutils.core import Command

# bootstrap setuptools for users that don't already have it
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

# from fipy.tests.testRunner import test, unittest

from setuptools.command.test import test as _test

def _TestClass(base):
    class _test(base):
        description = base.description + ", for FiPy and its examples"

        # List of option tuples: long name, short name (None if no short
        # name), and help string.
        user_options = base.user_options + [
            ('inline', None, "run FiPy with inline compilation enabled"),
            ('pythoncompiled=', None, "directory in which to put weave's work product"),
            ('Trilinos', None, "run FiPy using Trilinos solvers"),
            ('Pysparse', None, "run FiPy using Pysparse solvers (default)"),
            ('all', None, "run all non-interactive FiPy tests (default)"),
            ('really-all', None, "run *all* FiPy tests (including those requiring user input)"),
            ('examples', None, "test FiPy examples"),
            ('modules', None, "test FiPy code modules"),
            ('viewers', None, "test FiPy viewer modules (requires user input)"),
            ('cache', None, "run FiPy with Variable caching"),
            ('no-cache', None, "run FiPy without Variable caching"),
           ]


        def initialize_options(self):
            base.initialize_options(self)
            
            self.all = False
            self.really_all = False
            self.examples = False
            self.modules = False
            self.viewers = False
            
            self.inline = False
            self.pythoncompiled = None
            self.cache = False
            self.no_cache = True
            self.Trilinos = False
            self.Pysparse = False

        def finalize_options(self):
            noSuiteOrModule = (self.test_suite is None 
                               and self.test_module is None)
                
            base.finalize_options(self)
            
            if noSuiteOrModule:
                self.test_args.remove(self.distribution.test_suite)
                
            if not (self.examples or self.modules or self.viewers):
                self.all = True
            if self.all or self.really_all:
                self.examples = True
                self.modules = True
            if self.really_all:
                self.viewers = True
            
                
            if self.viewers:
                print "*" * 60
                print "*" + "".center(58) + "*"
                print "*" + "ATTENTION".center(58) + "*"
                print "*" + "".center(58) + "*"
                print "*" + "Some of the following tests require user interaction".center(58) + "*"
                print "*" + "".center(58) + "*"
                print "*" * 60
                
                self.test_args.append("fipy.viewers.testinteractive._suite")

            if self.modules:
                self.test_args.append("fipy.test._suite")
            
            if self.examples:
                self.test_args.append("examples.test._suite")

            if self.test_args and noSuiteOrModule:
                self.test_suite = "dummy"
                
        def run_tests(self):
            import sys
            if self.Trilinos:
                try:
                    ## The import scipy statement is added to allow
                    ## the --Trilinos tests to run without throwing a
                    ## segmentation fault. This is caused by weird
                    ## behavior in scipy and PyTrilinos depending on
                    ## the order in which modules are imported
                    try:
                        import scipy
                    except:
                        pass
                    import PyTrilinos
                except ImportError, a:
                    print >>sys.stderr, "!!! Trilinos library is not installed"
                    return

            if self.inline:
                try:
                    from scipy import weave
                except ImportError, a:
                    print >>sys.stderr, "!!! weave library is not installed"
                    return
                    
            if self.pythoncompiled is not None:
                import os
                os.environ['PYTHONCOMPILED'] = self.pythoncompiled

            base.run_tests(self)

    return _test                    
            
test = _TestClass(_test)

try:
    # we only need "unittest" if bitten is installed 
    # (and we're running as a bitten.slave)
    from bitten.util.testrunner import unittest as _unittest
    unittest = _TestClass(_unittest)
except ImportError, e:
    unittest = test



class build_docs (Command):

    description = "build the FiPy api documentation"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [('latex', None, "compile the LaTeX variant of the apis"),
                    ('html', None, "compile the HTML variant of the apis"),
                    ('guide', None, "compile the user guide"),
                    ('apis', None, "compile the programmer's reference"),
                    ('manual', None, "compile the manual"),
                    ('all', None, "compile both the LaTeX and HTML variants of the apis"),
                    ('webpage', None, "compile the html for the web page"),
                    ('upload', None, "upload webpages, documentation, and distributions to CTCMS website"),
                    ('uploadwww', None, "upload webpages to CTCMS website"),
                   ]


    def initialize_options (self):
        self.latex = 0
        self.html = 0
        self.guide = 0
        self.apis = 0
        self.manual = 0
        self.all = 0
        self.webpage = 0
        self.upload = 0
        self.uploadwww = 0
    # initialize_options()


    def finalize_options (self):
        if self.all:
            self.latex = 1
            self.manual = 1
            self.webpage = 1
            
        if self.manual:
            self.guide = 1
            self.apis = 1
            
    # finalize_options()

    def _initializeDirectory(self, dir, type = 'latex'):
        dir = os.path.join(dir, type)
        
        try:
            for root, dirs, files in os.walk(dir, topdown=False): 
                for name in files: 
                    os.remove(os.path.join(root, name)) 
                for name in dirs: 
                    os.rmdir(os.path.join(root, name)) 
            os.rmdir(dir)
        except:
            pass
            
        os.makedirs(dir)
        
    def _epydocFiles(self, module, dir = None, type = 'pdflatex'):
        dir = os.path.join(dir, type)
        
        import epydoc.cli
        epydoc.cli.cli(["--%s" % type, "--output", dir, 
                        "--no-private", "--show-imports", 
                        "--name", "FiPy", 
                        module])

    def _buildTeXAPIs(self):
        dir = os.path.join('documentation', 'manual', 'api')
        self._initializeDirectory(dir = dir, type = 'latex')
        dir = os.path.join(dir, 'latex')
        
        import epydoc.cli
        epydoc.cli.cli(["--latex", "--output", dir, 
                        "--graph=classtree", "--inheritance=listed", 
                        "--no-private", "--show-imports", 
                        "fipy/"])
        
        savedir = os.getcwd()
        try:
            
            os.chdir(os.path.join('documentation','manual', 'api','latex'))
            old = open('api.tex', 'r')
            new = open('api-rev.tex', 'w')
            
            new.write("""% This file is created automatically by:
% 	python setup.py build_docs --latex

""")
            
            import re
            
            mainModule = re.compile(r"^\\include{((fipy\.[^.-]*)-module)}")
            subModule = re.compile(r"^\\include{((fipy(\.[^.-]*)+)-module)}")

            for line in old:
                mainMatch = mainModule.match(line)
                subMatch = subModule.match(line)

                def stringInModule(s, name):
                    module = open(name, 'r')
                    functionLine = re.compile(s)
                    flag = False
                    for l in module:
                        if functionLine.search(l):
                            flag = True
                            break
                            
                    module.close()
                    
                    return flag


                if mainMatch:
                    moduleName = mainMatch.group(1) + ".tex"
                    if (stringInModule(r"\\section{Package", moduleName)
                        and not stringInModule("no chapter heading", moduleName)):
                        module = open(moduleName, 'r')
                        lines = []
                        
                        for l in module:
                            l = re.sub(r'\\section', r'\\chapter', l)
                            l = re.sub(r'\\subsection', r'\\section', l)
                            l = re.sub(r'\\subsubsection', r'\\subsection', l)
                            lines.append(l)
                            
                        module.close()
                        module = open(moduleName, 'w')
                        module.writelines(lines)
                        module.close()

                    if not stringInModule(r"\\section{(Functions|Variables|Class)", moduleName):
                        new.write("\\input{%s}\n" % mainMatch.group(1))
                    else:
                        new.write(line)
                elif subMatch:
                    ## epydoc tends to prattle on and on with empty module pages, so 
                    ## we eliminate all but those that actually contain something relevant.
                    moduleName = subMatch.group(1) + ".tex"
                    if not stringInModule(r"\\(sub)?section{(Functions|Variables|Class)", moduleName):
                        new.write("\\input{%s}\n" % subMatch.group(1))
                    else:
                        new.write(line)

            new.close()
            old.close()
        except Exception, e:
            print e
        
        os.chdir(savedir)
        
    def _translateTextFiles(self, source_dir = '.', destination_dir = '.', files = [], writer = None, settings = {}, ext = '.tex'):
        from docutils import core

        for file in files:

            destination_path = os.path.join(destination_dir, string.lower(file) + ext)
            try:
                os.makedirs(os.path.dirname(destination_path))
            except:
                pass
            source_path = os.path.join(source_dir, file + '.txt')

            core.publish_file(source_path= source_path,
                              destination_path = destination_path,
                              reader_name = 'standalone',
                              parser_name = 'restructuredtext',
                              writer = writer,
                              settings_overrides = settings)
                              
            translated = open(destination_path, 'r')
            lines = []
            
            for line in translated:
                import re
                line = re.sub(r'\\tableofcontents', r'% this automatically generated line conflicts with minitoc\r% \\tableofcontents', line)
                lines.append(line)
                
            translated.close()
            translated = open(destination_path, 'w')
            translated.writelines(lines)
            translated.close()

            # mark modification time of output file as mod time of reST file
            os.utime(destination_path, (os.path.getatime(source_path), os.path.getmtime(source_path)))

    def run (self):
        f = open(os.path.join('documentation','VERSION.txt'), 'w')
        f.write('.. |VERSION| replace:: ' + self.distribution.metadata.get_version())
        f.close()

        mainRestructuredTextFiles = {'article': 
                                         ['INSTALLATION',
                                          'README',
                                          'LICENSE',
                                          'DISCLAIMER',
                                          'examples/README'], 
                                     'startlower': 
                                         ['WINDOWS-INSTALLATION',
                                          'MACOSX-INSTALLATION',
                                          'examples/levelSet/electroChem/README']}
                                         
                                     
        secondaryRestructuredTextFiles = {'article': 
                                              ['CREDITS',
                                               'TALKS',
                                               'TODOLIST',
                                               'SVN',
                                               'EFFICIENCY',
                                               'VKML'],
                                          'startlower': 
                                              ['MAIL']}

        if self.latex:
            if self.apis:
                self._buildTeXAPIs()
                
            if self.guide:
                dir = os.path.join('documentation', 'manual', 'examples')
                self._initializeDirectory(dir = dir, type = 'latex')
                dir = os.path.join(dir, 'latex')
                               
                import epydoc.cli
                epydoc.cli.cli(["--latex", "--output", dir, 
                                "--no-private", "--show-imports", 
                                "examples/"])

        if self.html:
            dir = os.path.join('documentation', 'manual', 'api')
            self._initializeDirectory(dir = dir, type = 'html')
            self._epydocFiles(module = 'fipy/', dir = dir, type = 'html')

        if self.apis:
            # build the package/module/class example documentation
            
            dir = os.path.join('documentation', 'manual', 'tutorial')
            self._initializeDirectory(dir = dir, type = 'latex')
            dir = os.path.join(dir, 'latex')
##             self._initializeDirectory(dir = dir, type = 'pdflatex')
##             dir = os.path.join(dir, 'pdflatex')

            # to avoid a collision between the real fipy namespace
            # and the fictional fipy namespace we use for the illustration
            # we build the example documentation in a sub-process
            # from which we delete the "real" fipy namespace
            
            exec("""
import sys
if sys.modules.has_key('fipy'):
    del sys.modules['fipy']
    
if sys.modules.has_key('epydoc.uid'):
    sys.modules['epydoc.uid']._object_uids = {}
    sys.modules['epydoc.uid']._variable_uids = {}
    sys.modules['epydoc.uid']._name_to_uid = {}

import epydoc.cli
epydoc.cli.cli(["--latex", "--output", dir, 
                "--graph=classtree", 
                "--no-private", "--show-imports", 
                "documentation/manual/tutorial/fipy/"])
""")

        if self.guide or self.apis:
            savedir = os.getcwd()
            
            os.chdir(os.path.join('documentation','manual'))
                
            f = open('version.tex', 'w')
            f.write("% This file is created automatically by:\n")
            f.write("% 	python setup.py build_docs --manual\n\n")
            f.write("\\newcommand{\\Version}{" + self.distribution.metadata.get_version() + "}\n")
            f.close()
            
            from utils.includedLaTeXWriter import IncludedLaTeXWriter
            
            for key in mainRestructuredTextFiles.keys():
                self._translateTextFiles(files = mainRestructuredTextFiles[key],
                                         source_dir = '../..',
                                         writer = IncludedLaTeXWriter(),
                                         settings ={'use_latex_toc': True,
                                                    'footnote_references': 'superscript',
                                                    'table_style': 'nolines',
                                                    'documentclass': key,
                                                    'reference_label': "ref*"})

            for key in mainRestructuredTextFiles.keys():
                self._translateTextFiles(files = secondaryRestructuredTextFiles[key],
                                         source_dir = '..',
                                         writer = IncludedLaTeXWriter(),
                                         settings ={'use_latex_toc': True,
                                                    'footnote_references': 'superscript',
                                                    'table_style': 'booktabs',
                                                    'documentclass': key,
                                                    'reference_label': "ref*"})

            if self.guide:
                os.system("pdflatex fipy")
                os.system("bibtex fipy")
                os.system("makeindex fipy")
                os.system("pdflatex fipy")
                os.system("pdflatex fipy")
                
            if self.apis:
                os.system("pdflatex reference")
                os.system("bibtex reference")
                os.system("makeindex reference")
                os.system("pdflatex reference")
                os.system("pdflatex reference")

            os.chdir(savedir)

        if self.webpage:
            import tempfile
            tmp = tempfile.mkdtemp()
            dir = os.path.join('documentation', 'www')

            from utils.includedHTMLWriter import IncludedHTMLWriter
            
            print "main files"
            for key in mainRestructuredTextFiles.keys():
                self._translateTextFiles(files = mainRestructuredTextFiles[key],
                                         destination_dir = tmp,
                                         writer = IncludedHTMLWriter(),
                                         settings = {'initial_header_level' : 3,
                                                     'xml_declaration' : 0},
                                         ext = '.html')

            print "secondary files"
            for key in secondaryRestructuredTextFiles.keys():
                self._translateTextFiles(files = secondaryRestructuredTextFiles[key],
                                         source_dir = "documentation",
                                         destination_dir = tmp,
                                         writer = IncludedHTMLWriter(),
                                         settings = {'initial_header_level' : 3,
                                                     'xml_declaration' : 0},
                                         ext = '.html')

            import shutil
            for f in ['menu.html', 'meta.html', 'logo.html', 'extra.html']:
                shutil.copyfile(os.path.join(dir, f), os.path.join(tmp, f))
            shutil.move(os.path.join(tmp, 'readme.html'), os.path.join(tmp, 'index.html'))
            shutil.move(os.path.join(tmp, 'examples', 'levelSet', 'electroChem', 'readme.html'), os.path.join(tmp, 'electrochem.html'))
            
            print "merging files"
            os.system("/Library/WebServer/Documents/CSS/ctcmsWeb.py %s %s" % (tmp, dir))
            
            print "removing directories"
            for root, dirs, files in os.walk(tmp, topdown=False): 
                for name in files: 
                    os.remove(os.path.join(root, name)) 
                for name in dirs: 
                    os.rmdir(os.path.join(root, name)) 

        if self.uploadwww:
                 
            print "setting group and ownership of web pages..."
            os.system('chmod -R g+w documentation/www/')
            
            print "uploading web pages..."
            # The -t flag (implicit in -a) is suddenly causing problems
            # os.system('rsync -aLC -e ssh %s %s'%('documentation/www/', os.environ['FIPY_WWWHOST']))
            os.system('rsync -rlpgoDLC -e ssh %s %s'%('documentation/www/', os.environ['FIPY_WWWHOST']))

            print "activating web pages..."
            os.system(os.environ['FIPY_WWWACTIVATE'])
            
        if self.upload:
            print "setting permissions of manuals..."
            os.system('chmod -R g+w documentation/manual/fipy.pdf')
            os.system('chmod -R g+w documentation/manual/reference.pdf')
            
            print "linking manuals to `dist/`..."
            os.system('mkdir dist/')
            os.system('ln -f documentation/manual/fipy.pdf dist/fipy-%s.pdf'%self.distribution.metadata.get_version())
            os.system('ln -f documentation/manual/reference.pdf dist/reference-%s.pdf'%self.distribution.metadata.get_version())
            
            for name in ('.tar.gz', '.win32.zip'):
                file = 'dist/FiPy-%s%s'%(self.distribution.metadata.get_version(), name)
                print "setting permissions for %s ..."%file
                os.system('chmod -R g+w %s'%file)

            print "build products in `dist/` must be manually uploaded to MatForge"
            import webbrowser
            webbrowser.open("http://matforge.org/fipy/admin/general/downloader", autoraise=False)
            
            print "please update the current links, as appropriate"
            webbrowser.open("http://matforge.org/fipy/wiki/FiPyDownloadCurrent?action=edit", autoraise=False)
            webbrowser.open("http://matforge.org/fipy/wiki/FiPyManual?action=edit", autoraise=False)
            webbrowser.open("http://matforge.org/fipy/wiki/FiPyReference?action=edit", autoraise=False)
            


                
    # run()

class copy_script(Command):
    description = "copy an example script into a new editable file"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [
        # Select installation scheme and set base director(y|ies)
        ('From=', None,
         "path and file name containing script to copy"),
        ('To=', None,
         "path and file name to save script to")
     ]

    def initialize_options(self):
        self.From = None
        self.To = None

    def finalize_options(self):
        if self.From == None:
            raise "Please specify a '--From' input script file"
         
        if self.To == None:
            raise "Please specify a '--To' output script file"
            
        if os.path.exists(os.path.expanduser(self.To)):
            ans = "junk"
            
            while (len(ans) > 0) and ("yes".find(ans.lower()) is not 0) and ("no".find(ans.lower()) is not 0):
                ans = raw_input("The file '%s' already exists. Overwrite? [n] "%self.To)
                
            if ans is '':
                ans = 'no'
                
            if ("no".find(ans.lower()) is 0):
                self.To = raw_input("Please give a name for the ouput file: ")
                self.finalize_options()

    def run(self):
        import imp
        import fipy.tests.doctestPlus
        
        mod = imp.load_source("copy_script_module", self.From)
        script = fipy.tests.doctestPlus._getScript(name = "copy_script_module")
        
        script = "#!/usr/bin/env python\n\n## This script was derived from\n## '%s'\n\n%s"%(self.From, script)
        
        f = file(self.To, "w")
        f.write(script)
        f.close
        
        print "Script code exported from '%s' to '%s'"%(self.From, self.To)

class efficiency_test(Command):
    description = "run FiPy efficiency tests"
    
    user_options = [ ('minimumelements=', None, 'minimum number of elements'),
                     ('factor=', None, 'factor by which the number of elements is increased'),
                     ('inline', None, 'turn on inlining for the efficiency tests'),
                     ('cache', None, 'turn on variable caching'),
                     ('maximumelements=', None, 'maximum number of elements'),
                     ('sampleTime=', None, 'sampling interval for memory high-water'),
                     ('path=', None, 'directory to place output results in')]
    
    def initialize_options(self):
        self.factor = 10
        self.inline = 0
        self.cache = 0
        self.maximumelements = 10000
        self.minimumelements = 100
        self.sampleTime = 1
        self.path = None
        self.cases = ['examples/benchmarking/cahnHilliard.py', 'examples/benchmarking/superfill.py', 'examples/benchmarking/phaseImpingement.py', 'examples/benchmarking/mesh.py']
        
    def finalize_options(self):
        self.factor = int(self.factor)
        self.maximumelements = int(self.maximumelements)
        self.minimumelements = int(self.minimumelements)
        self.sampleTime = float(self.sampleTime)

    def run(self):

        import time
        import os
        
        for case in self.cases:
            print "case: %s" % case
            
            if self.path is None:
                testPath = os.path.split(case)[0]
            else:
                testPath = self.path
                
            if not os.access(testPath, os.F_OK):
                os.makedirs(testPath)
                
            testPath = os.path.join(testPath, '%s.dat' % os.path.split(case)[1])
            
            if not os.path.isfile(testPath):
                f = open(testPath, 'w')

                f.write("\t".join(["--inline", "--cache", "Date", "Elements", \
                                  "mesh (s)", "variables (s)", "terms (s)", \
                                  "solver (s)", "BCs (s)", "solve (s)", \
                                  "total (s)", "per step (s)", \
                                  \
                                  "mesh (KiB)", "variables (KiB)", \
                                  "terms (KiB)", "solver (KiB)", "BCs (KiB)", \
                                  "solve (KiB)", "max (KiB)", \
                                  "per element (KiB)"]))
                f.write("\n")
                f.flush()
            else:
                f = open(testPath, 'a')
            
            numberOfElements = self.minimumelements

            while numberOfElements <= self.maximumelements:
                print "\tnumberOfElements: %i" % numberOfElements
                
                cmd = [case, '--numberOfElements=%i' % numberOfElements]
                
                if self.inline:
                    cmd += ['--inline']
                    
                if self.cache:
                    cmd += ['--cache']
                else:
                    cmd += ['--no-cache']

                output = "\t".join([str(self.inline), str(self.cache), time.ctime(), str(numberOfElements)])
                
                timeCmd = cmd + ['--measureTime']
                w, r = os.popen4(' '.join(timeCmd))
                output += '\t' + ''.join(r.readlines()).strip()
                r.close()
                w.close()

                memCmd = cmd + ['--measureMemory', '--sampleTime=%f' % self.sampleTime]
                w, r = os.popen4(' '.join(memCmd))
                output += '\t' + ''.join(r.readlines()).strip()
                r.close()
                w.close()
                    
                f.write(output + '\n')
                f.flush()

                numberOfElements *= self.factor

            f.close()

try:            
    f = open('README.txt', 'r')
    long_description = '\n' + f.read() + '\n'
    f.close()
except IOError, e:
    long_description = ''
        
try:
    f = open('LICENSE.txt', 'r') 
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

dist = setup(	name = "FiPy",
        version = "2.0.2", 
        author = "Jonathan Guyer, Daniel Wheeler, & Jim Warren",
        author_email = "fipy@nist.gov",
        url = "http://www.ctcms.nist.gov/fipy/",
        license = license,
        description = "A finite volume PDE solver in Python",
        long_description = long_description,
        cmdclass = {
            'build_docs':build_docs,
            'test':test,
            'unittest':unittest,
            'copy_script': copy_script,
            'efficiency_test': efficiency_test
        },
        test_suite="fipy.test._suite",
        packages = find_packages(exclude=["examples", "examples.*", "utils", "utils.*"]),
        entry_points="""
            [fipy.viewers]
            gist = fipy.viewers.gistViewer:GistViewer
            gnuplot = fipy.viewers.gnuplotViewer:GnuplotViewer
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
    
    for pkg in ['scipy', 'matplotlib', 'gist', 'mayavi']:
        try:
            __import__(pkg)
        except ImportError, exc:
            opt.append(pkg)
        
    if len(opt) > 0:
        print "------"
        print "The optional module(s) " + str(opt) + " cannot be loaded."
        print "FiPy will have improved capabilities if these modules are installed."
