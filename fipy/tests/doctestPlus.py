#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "doctestPlus.py"
 #                                    created: 10/27/04 {9:14:53 AM} 
 #                                last update: 12/9/04 {8:28:45 PM} 
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
 # protection and is in the public domain.  doctest.py <2>
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
 #  2004-10-27 JEG 1.0 original
 # ###################################################################
 ##

import sys
import doctest

from lateImportTest import LateImportTestCase, LateImportTestSuite

def getScript(name = '__main__'):
    return doctest.testsource(sys.modules.get(name), "")
    
class LateImportDocTestCase(LateImportTestCase):
    def getTestSuite(self, module):
        return doctest.DocTestSuite(module)    

class LateImportDocTestSuite(LateImportTestSuite):
    def __init__(self, testModuleNames = (), docTestModuleNames = (), base = '__main__'):
        LateImportTestSuite.__init__(self, testModuleNames = testModuleNames, base = base)
        self.addDocTestModules(moduleNames = docTestModuleNames, base = base)
    
    def addDocTestModules(self, moduleNames = (), base = '__main__'):
        for moduleName in moduleNames:
            self.addTestModule(moduleName = moduleName, base = base, testClass = LateImportDocTestCase)

