#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 9/3/04 {10:40:08 PM} 
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

import doctest
import fipy.tests.testProgram
import unittest

import advection.advectionTerm
import advection.higherOrderAdvectionTerm
import distanceFunction.distanceVariable
import surfactant.surfactantVariable
import surfactant.adsorbingSurfactantEquation
import surfactant.convectionCoeff
import distanceFunction.levelSetDiffusionVariable
import electroChem.test

def suite():
    theSuite = unittest.TestSuite()

    theSuite.addTest(doctest.DocTestSuite(advection.advectionTerm))
    theSuite.addTest(doctest.DocTestSuite(advection.higherOrderAdvectionTerm))
    theSuite.addTest(doctest.DocTestSuite(surfactant.surfactantVariable))
    theSuite.addTest(doctest.DocTestSuite(surfactant.convectionCoeff))
    theSuite.addTest(doctest.DocTestSuite(surfactant.adsorbingSurfactantEquation))
    theSuite.addTest(doctest.DocTestSuite(distanceFunction.distanceVariable))
    theSuite.addTest(doctest.DocTestSuite(distanceFunction.levelSetDiffusionVariable))
    theSuite.addTest(electroChem.test.suite())
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')
