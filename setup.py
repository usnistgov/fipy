#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "setup.py"
 #                                    created: 4/6/04 {1:24:29 PM} 
 #                                last update: 11/2/04 {2:44:01 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James Warren <jwarren@nist.gov>
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
import string

from distutils.core import setup
from distutils.core import Command

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
                    ('upload', None, "upload webpages to lurch:/u/WWW/wd15/fipy")
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
    # initialize_options()


    def finalize_options (self):
	if self.all:
	    self.latex = 1
	    self.html = 1
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
	
    def _epydocFiles(self, module, dir = None, type = 'latex'):
	dir = os.path.join(dir, type)
	
        command = "epydoc --" + type + " --output " + dir + " --name FiPy " + module
        
	os.system(command)

    def _buildTeXAPIs(self):
	dir = os.path.join('documentation', 'manual', 'api')
	self._initializeDirectory(dir = dir, type = 'latex')
	self._epydocFiles(module = 'fipy/', dir = dir, type = 'latex')
	
        savedir = os.getcwd()
        try:
            
            os.chdir(os.path.join('documentation','manual'))
            f = open('api.tex', 'w')
            f.write("% This file is created automatically by:\n")
            f.write("% 	python setup.py build_doc --latex\n\n")
            for root, dirs, files in os.walk(os.path.join('api','latex'), topdown=True):
                
                if 'api.tex' in files:
                    files.remove('api.tex')
		    
		if 'fipy-module.tex' in files:
		    files.remove('fipy-module.tex')

                
                ## Added because linux does not sort files in the same order
                files.sort()
                
##                 for module in modules[::-1]:
##                     formattedModule = string.replace(module,'/','.') + '-module.tex'
##                     if formattedModule in files:
##                         files.remove(formattedModule)
##                         files.insert(0, formattedModule)

		import re
		mainModule = re.compile(r"(fipy\.[^.-]*)-module\.tex")
		subModule = re.compile(r"(fipy(\.[^.-]*)+)-module\.tex")
                for name in files:
		    mainMatch = mainModule.match(name)
		    if mainMatch:
			f.write("\\chapter{Module " + mainMatch.group(1) + "}\n")
			
		    subMatch = subModule.match(name)
		    if subMatch:
			module = open(os.path.join(root, name))
			
			# epydoc tends to prattle on and on with empty module pages, so 
			# we eliminate all but those that actually contain something relevant.
			functionLine = re.compile(r"\\subsection{(Functions|Variables|Class)")
			keepIt = False
			for line in module:
			    if functionLine.search(line):
				keepIt = True
				break
				
			module.close
			if not keepIt:
			    continue
			
		    split = os.path.splitext(name)
		    if split[1] == ".tex":
			f.write("\\input{" + os.path.join(root, os.path.splitext(name)[0]) + "}\n\\newpage\n")

            f.close()
        except:
            pass
        
        os.chdir(savedir)
	
    def _translateTextFiles(self, source_dir = '.', destination_dir = '.', files = [], writer = None, settings = {}, ext = '.tex'):
	from docutils import core

        for file in files:

            destination_path = os.path.join(destination_dir, string.lower(file) + ext)
            source_path = os.path.join(source_dir, file + '.txt')

            core.publish_file(source_path= source_path,
                              destination_path = destination_path,
                              reader_name = 'standalone',
                              parser_name = 'restructuredtext',
                              writer = writer,
                              settings_overrides = settings)

    def run (self):

        restructuredTextFiles = ['INSTALLATION',
                                 'README',
                                 'CREDITS',
                                 'TALKS',
                                 'TODOLIST',
				 'LICENSE',
				 'DISCLAIMER',
				 'MAIL',
				 'examples/README']
        
	if self.latex:
	    if self.apis:
		self._buildTeXAPIs()
		
	    if self.guide:
		dir = os.path.join('documentation', 'manual', 'examples')
		self._initializeDirectory(dir = dir, type = 'latex')
		for module in ['examples/diffusion/',
			       'examples/convection/',
			       'examples/phase/',
			       'examples/levelSet/',
			       'examples/elphf/',
			       'examples/cahnHilliard/'
			       ]:
		    self._epydocFiles(module = module, dir = dir, type = 'latex')


	if self.html:
	    dir = os.path.join('documentation', 'manual', 'api')
	    self._initializeDirectory(dir = dir, type = 'html')
	    self._epydocFiles(module = 'fipy/', dir = dir, type = 'html')

	if self.guide or self.apis:
	    savedir = os.getcwd()
	    
##	    try:
            os.chdir(os.path.join('documentation','manual'))
		
            f = open('version.tex', 'w')
            f.write("% This file is created automatically by:\n")
            f.write("% 	python setup.py build_doc --manual\n\n")
            f.write("\\newcommand{\\Version}{" + self.distribution.metadata.get_version() + "}\n")
            f.close()
            
	    from utils.includedLaTeXWriter import IncludedLaTeXWriter
	    
            self._translateTextFiles(files = restructuredTextFiles,
                                     source_dir = '../..',
                                     writer = IncludedLaTeXWriter(),
                                     settings ={'use_latex_toc': True,
                                                'footnote_references': 'superscript'})


	    if self.guide:
		os.system("pdflatex fipy")
                os.system("bibtex fipy")
		os.system("makeindex fipy")
		os.system("pdflatex fipy")
		os.system("pdflatex fipy")
		
	    if self.apis:
		os.system("pdflatex reference")
		os.system("makeindex reference")
		os.system("pdflatex reference")

##	    except:
##		pass
	    os.chdir(savedir)

        if self.webpage:
            dir = os.path.join('documentation', 'www')

	    from utils.includedHTMLWriter import IncludedHTMLWriter
	    
            self._translateTextFiles(files = restructuredTextFiles,
                                     destination_dir = dir,
                                     writer = IncludedHTMLWriter(),
                                     settings = {'initial_header_level' : 3,
                                                 'stylesheet' : 'ctcms.css',
						 'xml_declaration' : 0},
                                     ext = '.html')

            import shutil
            shutil.move(os.path.join(dir, 'readme.html'), os.path.join(dir, 'index.html'))

        if self.upload:

	    wwwhost = 'dromio.nist.gov:/u/WWW/wd15/fipy/'

	    print "setting group and ownership of manuals..."
	    os.system('chgrp -R pfm documentation/manual/fipy.pdf')
	    os.system('chmod -R g+w documentation/manual/reference.pdf')
	    os.system('chmod -R g+w documentation/manual/reference.pdf')
	    os.system('chgrp -R pfm documentation/manual/fipy.pdf')
	    os.system('chmod -R g+w dist/FiPy-%s.tar.gz'%self.distribution.metadata.get_version())
	    os.system('chgrp -R pfm dist/FiPy-%s.tar.gz'%self.distribution.metadata.get_version())
	    
	    print "linking manuals to website..."
	    os.system('ln -sf ../../manual/fipy.pdf documentation/www/download/')
	    os.system('ln -sf ../../manual/reference.pdf documentation/www/download/')
	    os.system('ln -sf ../../../dist/FiPy-%s.tar.gz documentation/www/download/'%self.distribution.metadata.get_version())
	    
	    print "setting group and ownership of web pages..."
	    os.system('chgrp -R pfm documentation/www/')
	    os.system('chmod -R g+w documentation/www/')
	    
	    print "uploading web pages..."
	    os.system('rsync -aLC -e ssh %s %s'%('documentation/www/',wwwhost))

	    print "activating web pages..."
	    os.system('ssh dromio.nist.gov "ssh lurch updatewww"')

                
    # run()

class test(Command):
    description = "test FiPy and its examples"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [('inline', None, "run FiPy with inline compilation enabled"),
		    ('all', None, "run all FiPy tests (default)"),
		    ('examples', None, "test FiPy examples"),
		    ('modules', None, "test FiPy code modules"),
		    ('terse', None, "give limited output during tests"),
		   ]


    def initialize_options (self):
	self.inline = False
	self.verbosity = 0
	self.terse = False
	self.all = False
	self.doExamples = True
	self.doModules = True
	self.examples = False
	self.modules = False

    def finalize_options (self):
	if self.verbose:
	    self.verbosity = 2
	if self.terse:
	    self.verbosity = 1
	if self.all:
	    self.examples = True
	    self.modules = True
	if self.examples and not self.modules:
	    self.doModules = False
	if self.modules and not self.examples:
	    self.doExamples = False
	
    def run (self):
	import unittest
	theSuite = unittest.TestSuite()
	
	if self.doModules:
	    import fipy.test
	    theSuite.addTest(fipy.test.suite())
	
	if self.doExamples:
	    import examples.test
	    theSuite.addTest(examples.test.suite())
	
	testRunner = unittest.TextTestRunner(verbosity=self.verbosity)
	result = testRunner.run(theSuite)
	
	import sys
	sys.exit(not result.wasSuccessful())

	    
long_description = """
A finite volume PDE solver in Python.

The authors and maintainers of this package are:
    
Daniel Wheeler <daniel.wheeler@nist.gov>
Jonathan Guyer <guyer@nist.gov>
Jim Warren <jwarren@nist.gov>
"""

f = open('LICENSE.txt', 'r') 
license = '\n' + f.read() + '\n'
f.close()

dist = setup(	name = "FiPy",
	version = "0.1",
	author = "Jonathan Guyer, Daniel Wheeler, & Jim Warren",
	author_email = "guyer@nist.gov",
	url = "http://ctcms.nist.gov/fipy/",
	license = license,
	description = "A finite volume PDE solver in Python",
	long_description = long_description,
	cmdclass = {
	    'build_docs':build_docs,
	    'test':test
	},
	packages = ['fipy', 
			'fipy.boundaryConditions',
			'fipy.equations',
			'fipy.iterators',
			'fipy.meshes',
			    'fipy.meshes.common',
			    'fipy.meshes.numMesh',
			    'fipy.meshes.pyMesh',
			'fipy.models',
			    'fipy.models.cahnHilliard',
			    'fipy.models.elphf',
			    'fipy.models.levelSet',
				'fipy.models.levelSet.advection',
				'fipy.models.levelSet.distanceFunction',
				'fipy.models.levelSet.electroChem',
				'fipy.models.levelSet.surfactant',
			    'fipy.models.phase',
				'fipy.models.phase.phase',
				'fipy.models.phase.temperature',
				'fipy.models.phase.theta',
			'fipy.solvers',
			'fipy.terms',
			'fipy.tests',
			'fipy.tools',
			    'fipy.tools.dimensions',
			    'fipy.tools.inline',
			'fipy.variables',
			'fipy.viewers'
	]
)

if 'install' in dist.commands:
    req = []
    
    try:
	import Numeric
    except ImportError, exc:
	req.append('Numeric')
	
    try:
	import spmatrix
    except ImportError, exc:
	req.append('PySparse')
	
    if len(req) > 0:
	print "!!!!!!"
	print "The required module(s) " + str(req) + " cannot be loaded."
	print "FiPy will not work properly until these modules are installed."

    opt = []
    
    try:
	import weave
    except ImportError, exc:
	opt.append('weave')
	
    try:
	import gmsh
    except ImportError, exc:
	opt.append('gmsh')
	
    if len(opt) > 0:
	print "------"
	print "The optional module(s) " + str(opt) + " cannot be loaded."
	print "FiPy will have improved capabilities if these modules are installed."
	
    if not os.environ.has_key('FIPY_NOGIST'):
	try:
	    import gist
	except ImportError, exc:
	    print "------"
	    print "The Pygist package cannot be loaded."
	    print "If you choose not to install this package, be sure to create an " + \
	    "environment variable FIPY_NOGIST to allow the test suite to run to completion."
	
    if not os.environ.has_key('FIPY_NOPYX'):
	try:
	    import pyx
	except ImportError, exc:
	    print "------"
	    print "The PyX package cannot be loaded."
	    print "If you choose not to install this package, be sure to create an " + \
	    "environment variable FIPY_NOPYX to allow the test suite to run to completion."
	
