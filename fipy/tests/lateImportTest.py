## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "lateImportTest.py"
 #                                    created: 12/9/04 {8:26:13 PM} 
 #                                last update: 4/1/05 {2:48:40 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  lateImportTest.py
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
 #  2004-12-09 JEG 1.0 original
 # ###################################################################
 ##

import unittest
 
class LateImportTestCase(unittest.TestCase):
    def __str__(self):
        return "import %s" % self.moduleName

    def post__init__(self, moduleName, suite):
        self.moduleName = moduleName
        self.suite = suite
        
    def setUp(self):
        """
        See documentation of `__import__` for why
        this ugly hack is necessary
        """
        module = __import__(self.moduleName)
        components = self.moduleName.split('.')
        for component in components[1:]:
            module = getattr(module, component)
            
        self.suite.addTest(self.getTestSuite(module = module))
        
    def getTestSuite(self, module):
        return module._suite()
        
    def runTest(self):
        return
    
class LateImportTestSuite(unittest.TestSuite):
    def __init__(self, testModuleNames = (), base = '__main__'):
        unittest.TestSuite.__init__(self)
        self.addTestModules(moduleNames = testModuleNames, base = base)
    
    def addTestModules(self, moduleNames = (), base = '__main__'):
        for moduleName in moduleNames:
            self.addTestModule(moduleName = moduleName, base = base)

    def addTestModule(self, moduleName, base = '__main__', testClass = LateImportTestCase):
        if base == '__main__':
            base = []
        else:
            base = base.split('.')[:-1]

        test = testClass()
        test.post__init__(".".join(base + [moduleName]), self)
        self.addTest(test)

