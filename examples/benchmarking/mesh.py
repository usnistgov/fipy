#!/usr/bin/env python

## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "mesh.py"
 #                                     created: 1/18/06 {4:01:30 PM}
 #                                 last update: 1/18/06 {5:14:13 PM}
 # Author: Jonathan Guyer
 # E-mail: <guyer@nist.gov>
 # Author: Daniel Wheeler
 # E-mail: daniel.wheeler@nist.gov
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are 
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 # 
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 1940-01-18 JEG 1.0 original
 # 
 # ########################################################################
 ##

r"""
This example benchmarks the speed and memory usage of creating a mesh. Run:
    
    $ python setup.py efficiency_test
"""
__docformat__ = 'restructuredtext'

from fipy.tools.parser import parse

from benchmarker import Benchmarker
bench = Benchmarker()

numberOfElements = parse('--numberOfElements', action = 'store', type = 'int', default = 100)

bench.start()

from fipy.tools import numerix
nx = int(numerix.sqrt(numberOfElements))
ny = nx
dx = 1.
dy = 1.

from fipy.meshes.grid2D import Grid2D
mesh = Grid2D(nx = nx, ny = nx, dx = dx, dy = dy)

bench.stop('mesh')

print bench.report(numberOfElements=numberOfElements)
