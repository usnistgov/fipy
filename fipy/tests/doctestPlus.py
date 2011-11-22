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
import time
import doctest

__all__ = ["execButNoTest", "register_skipper", "report_doctest_skips"]

_DocTestTimes = []

from fipy.tests.lateImportTest import _LateImportTestCase, _LateImportTestSuite

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
    
_doctestSkippers = list()

def register_skipper(flag, test, why, skipWarning=True):
    global _doctestSkippers
    
    skipper = _DoctestSkipper(flag=doctest.register_optionflag(flag),
                              test=test,
                              why=why,
                              skipWarning=skipWarning)
    _doctestSkippers.append(skipper)

def report_doctest_skips():
    global _doctestSkippers
    
    skips = list()
    for skipper in _doctestSkippers:
        if skipper.skipWarning and skipper.skipped:
            skips.append("Skipped %d doctest examples because %s" 
                         % (len(skipper.skipped), skipper.why))
    if len(skips) > 0:
        print >>sys.stderr, "!" * 79
        print >>sys.stderr, "\n".join(skips)
        print >>sys.stderr, "!" * 79
    
class _DoctestSkipper:
    def __init__(self, flag, test, why, skipWarning):
        self.flag = flag
        self.why = why
        self.test = test
        self.skipWarning = skipWarning
        self.skipped = list()
        
    def skipTest(self):
        if not hasattr(self, "hasFeature"):
            self.hasFeature = self.test()
        return not self.hasFeature
        
def _checkForSciPy():
    hasSciPy = True
    try:
        import scipy
    except Exception:
        hasSciPy = False
    return hasSciPy
    
register_skipper(flag="SCIPY",
                 test=_checkForSciPy,
                 why="the `scipy` package cannot be imported")

class _SelectiveDocTestParser(doctest.DocTestParser):
    """ 
    Custom doctest parser that adds support for skipping test examples
    """ 
    def parse(self, string, name='<string>'): 
        pieces = doctest.DocTestParser.parse(self, string, name) 
        
        return [piece for piece in pieces if not self._skipExample(piece)]
        
    def _skipExample(self, piece):
        global _doctestSkippers
        
        skip = False
        
        if isinstance(piece, doctest.Example):
            for skipper in _doctestSkippers:
                if (piece.options.get(skipper.flag, False) and skipper.skipTest()):
                    skip = True
                    skipper.skipped.append(piece)
                    break
        
        return skip

    
class _LateImportDocTestCase(_LateImportTestCase):
    def _getTestSuite(self, module):
        return doctest.DocTestSuite(module, 
                                    test_finder=doctest.DocTestFinder(parser=_SelectiveDocTestParser()),
                                    setUp=self._setUp, tearDown=self._tearDown)
        
        
    @staticmethod
    def _setUp(docTestObj):
        docTestObj._startTime = time.time()
        docTestObj._endTime = docTestObj._startTime

    @staticmethod
    def _tearDown(docTestObj):
        docTestObj._endTime = time.time()
        _DocTestTimes.append((docTestObj._endTime - docTestObj._startTime, docTestObj.name))


class _LateImportDocTestSuite(_LateImportTestSuite):
    def __init__(self, testModuleNames=(), 
                 docTestModuleNames=(), 
                 base='__main__'):
        _LateImportTestSuite.__init__(self, testModuleNames = testModuleNames, base = base)
        self._addDocTestModules(moduleNames=docTestModuleNames, base=base)
    
    def _addDocTestModules(self, moduleNames=(), base='__main__'):
        for moduleName in moduleNames:
            self._addTestModule(moduleName=moduleName, base=base, testClass=_LateImportDocTestCase)

