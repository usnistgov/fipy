#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "setup.py"
 #                                    created: 4/6/04 {1:24:29 PM} 
 #                                last update: 6/7/04 {2:08:08 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from distutils.core import setup
import glob
import os
import string
from distutils.core import Command

class build_docs (Command):

    description = "build the FiPy api documentation"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [('latex', None, "compile the LaTeX variant of the apis"),
		    ('html', None, "compile the HTML variant of the apis"),
		    ('manual', None, "compile the manual"),
		    ('all', None, "compile both the LaTeX and HTML variants of the apis")
		   ]


    def initialize_options (self):
	self.latex = 0
	self.html = 0
	self.manual = 0
	self.all = 0

    # initialize_options()


    def finalize_options (self):
	if self.all:
	    self.latex = 1
	    self.html = 1
	    self.manual = 1
	    
    # finalize_options()

    def _buildTexFiles(self, modules, type = 'latex'):
	dir = os.path.join('documentation', 'manual', 'api', type)
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
        
        command = "epydoc --" + type + " --output " + dir + " --name FiPy --docformat restructuredtext "
        for module in modules:
            command = command + module + "/ "
        
	os.system(command)

    def _buildApiFile(self, modules):
        savedir = os.getcwd()
        try:
            
            os.chdir(os.path.join('documentation','manual'))
            f = open('api.tex', 'w')
            f.write("% This file is created automatically by:\n")
            f.write("% 	python setup.py build_doc --latex\n\n")
            for root, dirs, files in os.walk(os.path.join('api','latex'), topdown=True):
                
                if 'api.tex' in files:
                    files.remove('api.tex')
                
                ## Added because linux does not sort files in the same order
                files.sort()
                
                for module in modules[::-1]:
                    formattedModule = string.replace(module,'/','.') + '-module.tex'
                    if formattedModule in files:
                        files.remove(formattedModule)
                        files.insert(0, formattedModule)

                for name in files:
                    f.write("\\include{" + os.path.join(root, os.path.splitext(name)[0]) + "}\n")

            f.close()
        except:
            pass
        
        os.chdir(savedir)

    def run (self):
        modules = ['fipy',
                   'examples/diffusion',
                   'examples/convection',
                   'examples/phase/anisotropy',
                   'examples/phase/impingement',
                   'examples/phase/missOrientation']

	if self.latex:
             
	    self._buildTexFiles(modules = modules)
            self._buildApiFile(modules = modules)

	if self.html:
	    self._buildTexFiles(modules = modules, type = 'html')

	if self.manual:
	    savedir = os.getcwd()
	    try:
		os.chdir(os.path.join('documentation','manual'))
		os.system("pdflatex fipy.tex")
	    except:
		pass
	    os.chdir(savedir)

    # run()

# class x

long_description = """
A finite volume PDE solver in Python.

The authors and maintainers of this package are:
    
Daniel Wheeler <daniel.wheeler@nist.gov>
Jonathan Guyer <guyer@nist.gov>
Jim Warren <jwarren@nist.gov>
"""

setup(	name = "FiPy",
	version = "0.1",
	author = "Jonathan Guyer, Daniel Wheeler, & Jim Warren",
	author_email = "guyer@nist.gov",
	url = "http://ctcms.nist.gov",
	description = "A finite volume PDE solver in Python",
	long_description = long_description,
	
	cmdclass = {'build_docs':build_docs},
	packages = ['fipy', 
			'fipy.boundaryConditions',
			'fipy.equations',
			'fipy.iterators',
			'fipy.meshes',
			    'fipy.meshes.numMesh',
			    'fipy.meshes.pyMesh',
			'fipy.models',
			    'fipy.models.elphf',
			    'fipy.models.levelSet',
			    'fipy.models.phase',
				'fipy.models.phase.phase',
				'fipy.models.phase.temperature',
			'fipy.solvers',
			'fipy.terms',
			'fipy.tests',
			'fipy.tools',
			    'fipy.tools.dimensions',
			    'fipy.tools.inline',
			    'fipy.tools.profiler',
			'fipy.variables',
			'fipy.viewers'
	]
)
