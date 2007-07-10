#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "parallel_test.py"
 #                                    created: 4/6/04 {1:24:29 PM} 
 #                                last update: 2/12/07 {3:03:34 PM} 
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

# This file is necessary, because mpirun adds extraneous arguments to the 
# scripts it invokes, which breaks the option-processing of Distutils, which
# tries to handle all arguments and quits if it cannot. Thus, I haven't been
# able to get "mpirun -np NUM_PROC setup.py parallel_test --Trilinos" to work
# the way it should.

# should be invoked with mpirun -np NUM_PROC parallel_test.py --Trilinos

# Runs all non-interactive tests, with verbosity 2.

# All but one processor will fail the Trilinos matrix tests; the other tests
# should give the same results on all processors.

import unittest
theSuite = unittest.TestSuite()

import fipy.test
theSuite.addTest(fipy.test._suite())

import examples.test
theSuite.addTest(examples.test._suite())
    
testRunner = unittest.TextTestRunner(verbosity=2)
result = testRunner.run(theSuite)
    
import sys
sys.exit(not result.wasSuccessful())
