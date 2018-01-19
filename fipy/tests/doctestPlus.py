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
import inspect

__all__ = ["execButNoTest", "register_skipper", "report_skips", "testmod"]

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
    """Create a new doctest option flag for skipping tests

    Parameters
    ----------
    flag : str
      Name of the option flag
    test : function
      A function which should return `True` if the test should be run
    why : str
      Explanation for why the test was skipped (to be used in a string
      "``Skipped %%(count)d doctest examples because %%(why)s``")
    skipWarning : bool
      Whether or not to report on tests skipped by this flag (default `True`)
    """
    global _doctestSkippers

    skipper = _DoctestSkipper(flag=doctest.register_optionflag(flag),
                              test=test,
                              why=why,
                              skipWarning=skipWarning)
    _doctestSkippers.append(skipper)

def report_skips():
    """Print out how many doctest examples were skipped due to flags
    """
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

def testmod(m=None, name=None, globs=None, verbose=None,
            report=True, optionflags=0, extraglobs=None,
            raise_on_error=False, exclude_empty=False):
    """Test examples in the given module.  Return (#failures, #tests).

    Largely duplicated from :func:`doctest.testmod`, but using
    :class:`_SelectiveDocTestParser`.

    Test examples in docstrings in functions and classes reachable
    from module m (or the current module if m is not supplied), starting
    with m.__doc__.

    Also test examples reachable from dict m.__test__ if it exists and is
    not None.  m.__test__ maps names to functions, classes and strings;
    function and class docstrings are tested even if the name is private;
    strings are tested directly, as if they were docstrings.

    Return (#failures, #tests).

    See help(doctest) for an overview.

    Optional keyword arg "name" gives the name of the module; by default
    use m.__name__.

    Optional keyword arg "globs" gives a dict to be used as the globals
    when executing examples; by default, use m.__dict__.  A copy of this
    dict is actually used for each docstring, so that each docstring's
    examples start with a clean slate.

    Optional keyword arg "extraglobs" gives a dictionary that should be
    merged into the globals that are used to execute examples.  By
    default, no extra globals are used.  This is new in 2.4.

    Optional keyword arg "verbose" prints lots of stuff if true, prints
    only failures if false; by default, it's true iff "-v" is in sys.argv.

    Optional keyword arg "report" prints a summary at the end when true,
    else prints nothing at the end.  In verbose mode, the summary is
    detailed, else very brief (in fact, empty if all tests passed).

    Optional keyword arg "optionflags" or's together module constants,
    and defaults to 0.  This is new in 2.3.  Possible values (see the
    docs for details):

        DONT_ACCEPT_TRUE_FOR_1
        DONT_ACCEPT_BLANKLINE
        NORMALIZE_WHITESPACE
        ELLIPSIS
        SKIP
        IGNORE_EXCEPTION_DETAIL
        REPORT_UDIFF
        REPORT_CDIFF
        REPORT_NDIFF
        REPORT_ONLY_FIRST_FAILURE

    as well as FiPy's flags

        GMSH
        SCIPY
        TVTK
        SERIAL
        PARALLEL
        PROCESSOR_0
        PROCESSOR_0_OF_2
        PROCESSOR_1_OF_2
        PROCESSOR_0_OF_3
        PROCESSOR_1_OF_3
        PROCESSOR_2_OF_3

    Optional keyword arg "raise_on_error" raises an exception on the
    first unexpected exception or failure. This allows failures to be
    post-mortem debugged.
    """
    # If no module was given, then use __main__.
    if m is None:
        # DWA - m will still be None if this wasn't invoked from the command
        # line, in which case the following TypeError is about as good an error
        # as we should expect
        m = sys.modules.get('__main__')

    # Check that we were actually given a module.
    if not inspect.ismodule(m):
        raise TypeError("testmod: module required; %r" % (m,))

    # If no name was given, then use the module's name.
    if name is None:
        name = m.__name__

    # Find, parse, and run all tests in the given module.
    finder = doctest.DocTestFinder(exclude_empty=exclude_empty,
                                   parser=_SelectiveDocTestParser())

    if raise_on_error:
        runner = doctest.DebugRunner(verbose=verbose, optionflags=optionflags)
    else:
        runner = doctest.DocTestRunner(verbose=verbose, optionflags=optionflags)

    for test in finder.find(m, name, globs=globs, extraglobs=extraglobs):
        runner.run(test)

    if report:
        runner.summarize()
        report_skips()

    from fipy.tools import numerix
    printoptions = numerix.get_printoptions()
    if "legacy" in printoptions:
        numerix.set_printoptions(legacy="1.13")

    results = doctest.TestResults(runner.failures, runner.tries)

    if "legacy" in printoptions:
        numerix.set_printoptions(legacy=printoptions["legacy"])

    return results
