#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "test.py"
 #                                    created: 7/28/04 {11:18:34 AM} 
 #                                last update: 7/30/04 {7:16:33 PM} 
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
 # protection and is in the public domain.  test.py
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
 #  See the file "license.terms" for information on usage and  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-07-28 JEG 1.0 original
 # ###################################################################
 ##

def suite():
    import doctest
    import unittest
    
    theSuite = unittest.TestSuite()

    import input1D
    theSuite.addTest(doctest.DocTestSuite(input1D))
    import input1Ddimensional
    theSuite.addTest(doctest.DocTestSuite(input1Ddimensional))
    import input2D
    theSuite.addTest(doctest.DocTestSuite(input2D))
    import input2Dcorner
    theSuite.addTest(doctest.DocTestSuite(input2Dcorner))
    import input1Dphase
    theSuite.addTest(doctest.DocTestSuite(input1Dphase))
    import input1DphaseBinary
    theSuite.addTest(doctest.DocTestSuite(input1DphaseBinary))
    import input1DphaseQuaternary
    theSuite.addTest(doctest.DocTestSuite(input1DphaseQuaternary))
    import input1DphaseTernAndElectrons
    theSuite.addTest(doctest.DocTestSuite(input1DphaseTernAndElectrons))
    import input1DpoissonAllCharge
    theSuite.addTest(doctest.DocTestSuite(input1DpoissonAllCharge))
    import input1DpoissonLeftCharge
    theSuite.addTest(doctest.DocTestSuite(input1DpoissonLeftCharge))
    import input1DpoissonRightCharge
    theSuite.addTest(doctest.DocTestSuite(input1DpoissonRightCharge))

    return theSuite
    
if __name__ == '__main__':
    import fipy.tests.testProgram
    fipy.tests.testProgram.main(defaultTest='suite')
