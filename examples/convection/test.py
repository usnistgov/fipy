#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 12/7/04 {4:25:48 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

import unittest
import fipy.tests.testProgram

import doctest

import exponential1D.input
import exponential2D.input
import exponential1DBack.input
import powerLaw1D.input
import exponential1DSource.input
import exponential2D.tri2Dinput
import exponential1DSource.tri2Dinput

def suite():
    theSuite = unittest.TestSuite()

    theSuite.addTest(doctest.DocTestSuite(exponential1D.input))
    theSuite.addTest(doctest.DocTestSuite(exponential2D.input))
    theSuite.addTest(doctest.DocTestSuite(exponential1DBack.input))
    theSuite.addTest(doctest.DocTestSuite(powerLaw1D.input))
    theSuite.addTest(doctest.DocTestSuite(exponential1DSource.input))
    theSuite.addTest(doctest.DocTestSuite(exponential2D.tri2Dinput))
    theSuite.addTest(doctest.DocTestSuite(exponential1DSource.tri2Dinput))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')

            
            
