#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "doctestPlus.py"
 #
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
 # protection and is in the public domain.  doctestPlus.py
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
 # ###################################################################
 ##

import sys
import doctest

from lateImportTest import _LateImportTestCase, _LateImportTestSuite

def _getScript(name = '__main__'):
    module = sys.modules.get(name)
    # the syntax of doctest changed substantially between Python 2.3 and 2.4
    # <http://sourceforge.net/tracker/index.php?func=detail&aid=1120348&group_id=118428&atid=681141>
    if sys.version_info >= (2, 4):
        # Python 2.4 returns comments, too, and doesn't always end in a \n,
        # which chokes exec/compile. Arguably a bug in Python.
        # <http://sourceforge.net/tracker/index.php?func=detail&aid=1172785&group_id=5470&atid=105470>
        return doctest.testsource(module, module.__name__) + '\n'
    else:
        return doctest.testsource(module, "")
        
def execButNoTest(name='__main__'):
    module = sys.modules.get(name)
    
    # the syntax of doctest changed substantially between Python 2.3 and 2.4
    # <http://sourceforge.net/tracker/index.php?func=detail&aid=1120348&group_id=118428&atid=681141>
    if sys.version_info >= (2, 4):
        tests = doctest.DocTestFinder().find(module)
        tests = [doctest.script_from_examples(t.docstring) for t in tests]
        
        # Python 2.4 returns comments, too, and doesn't always end in a \n,
        # which chokes exec/compile. Arguably a bug in Python.
        # <http://sourceforge.net/tracker/index.php?func=detail&aid=1172785&group_id=5470&atid=105470>
        tests = [t + '\n' for t in tests]
    else:
        tests = [doc for (dummy, doc, dummy, dummy) in doctest._find_tests(module, "")]
        tests = [doctest._extract_examples(t) for t in tests]
	tests = ["\n".join([source for source, expect, dummy in t]) for t in tests]

    if not tests:
        raise ValueError("no tests found")

    for t in tests:
        exec t
        
class _LateImportDocTestCase(_LateImportTestCase):
    def _getTestSuite(self, module):
        return doctest.DocTestSuite(module)

class _LateImportDocTestSuite(_LateImportTestSuite):
    def __init__(self, testModuleNames=(), 
                 docTestModuleNames=(), 
                 base='__main__'):
        _LateImportTestSuite.__init__(self, testModuleNames = testModuleNames, base = base)
        self._addDocTestModules(moduleNames=docTestModuleNames, base=base)
    
    def _addDocTestModules(self, moduleNames=(), base='__main__'):
        for moduleName in moduleNames:
            self._addTestModule(moduleName=moduleName, base=base, testClass=_LateImportDocTestCase)

