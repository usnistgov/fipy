#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:00:39 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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

"""Test steady-state diffusion solutions
"""

from fipy.tests.testBase import TestBase
import Numeric



class TestLevelSetBase(TestBase):
##    def testinitialEvaluatedCellIDs(self):
##        self.assertArrayEqual(self._testInitialEvaluatedIDs, self.system.getInitialEvaluatedIDs())
    def getInitialEvaluatedValues(self):
        self.eqn.resetCells()
        self.eqn.setInitialEvaluatedCells()
        return self.eqn.getVar().getNumericValue()
    
    def getInitialTrialCellIDs(self):
        self.eqn.resetCells()
        self.eqn.setInitialEvaluatedCells()
        trialCellIDs = self.eqn.getInitialTrialCells()
        return trialCellIDs

    def getFinalValues(self):
        self.eqn.solve()
        return self.eqn.getVar().getNumericValue()

    def getInitialTrialValues(self):
        self.eqn.resetCells()
        self.eqn.setInitialEvaluatedCells()
        self.eqn.getInitialTrialCells()
        return self.eqn.getVar().getNumericValue()
    
    
    def testInitialEvaluatedValues(self):
        self.assertArrayWithinTolerance(self._testInitialEvaluatedValues, self.getInitialEvaluatedValues())

    def testInitialTrialCellIDs(self):
        self.assertArrayEqual(Numeric.sort(self._testInitialTrialCellIDs), Numeric.sort(self.getInitialTrialCellIDs())) 
        
    def testInitialTrialValues(self):
        self.assertArrayWithinTolerance(self._testInitialTrialValues, self.getInitialTrialValues())

    def testFinalValues(self):
        self.assertArrayWithinTolerance(self._testFinalValues, self.getFinalValues())
        
    def testResult(self):
        pass

    def evalCell(self, phix, phiy, dx, dy):
        aa = dy**2 + dx**2
        bb = -2 * ( phix * dy**2 + phiy * dx**2)
        cc = dy**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dy**2
        sqr = Numeric.sqrt(bb**2 - 4. * aa * cc)
        return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
        
        
