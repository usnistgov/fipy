#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "testClass.py"
 #
 #  Author: Jonathan Guyer   <guyer@nist.gov>
 #  Author: Daniel Wheeler   <daniel.wheeler@nist.gov>
 #  Author: James Warren     <jwarren@nist.gov>
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

def _TestClass(base):
    class _test(base):
        description = str(base.description) + ", for FiPy and its examples"

        # List of option tuples: long name, short name (None if no short
        # name), and help string.
        user_options = base.user_options + [
            ('inline', None, "run FiPy with inline compilation enabled"),
            ('pythoncompiled=', None, "directory in which to put weave's work product"),
            ('Trilinos', None, "run FiPy using Trilinos solvers"),
            ('Pysparse', None, "run FiPy using Pysparse solvers (default)"),
            ('trilinos', None, "run FiPy using Trilinos solvers"),
            ('pysparse', None, "run FiPy using Pysparse solvers (default)"),
            ('scipy', None, "run FiPy using SciPy solvers"),
            ('Scipy', None, "run FiPy using SciPy solvers"),
            ('no-pysparse',None, "run FiPy without using the Pysparse solvers"),
            ('pyamg',None, "run FiPy without using the PyAMG solvers"),
            ('all', None, "run all non-interactive FiPy tests (default)"),
            ('really-all', None, "run *all* FiPy tests (including those requiring user input)"),
            ('examples', None, "test FiPy examples"),
            ('modules', None, "test FiPy code modules"),
            ('viewers', None, "test FiPy viewer modules (requires user input)"),
            ('cache', None, "run FiPy with Variable caching"),
            ('no-cache', None, "run FiPy without Variable caching"),
            ('timetests=', None, "file in which to put time spent on each test"),
            ('skfmm', None, "run FiPy using the Scikit-fmm level set solver (default)"),
            ('lsmlib', None, "run FiPy using the LSMLIB level set solver (default)"),
           ]


        def initialize_options(self):
            base.initialize_options(self)

            self.all = False
            self.really_all = False
            self.examples = False
            self.modules = False
            self.viewers = False

            self.inline = False
            self.pythoncompiled = None
            self.cache = False
            self.no_cache = True
            self.Trilinos = False
            self.Pysparse = False
            self.trilinos = False
            self.pysparse = False
            self.no_pysparse = False
            self.pyamg = False
            self.scipy = False
            self.timetests = None
            self.skfmm = False
            self.lsmlib = False

        def finalize_options(self):
            noSuiteOrModule = (self.test_suite is None
                               and self.test_module is None)

            base.finalize_options(self)

            if noSuiteOrModule:
                # setuptools completely changed how it uses test_suite and test_args with v. 18.0
                # we do our best to keep it confused
                self.test_suite = None
                
            if not (self.examples or self.modules or self.viewers):
                self.all = True
            if self.all or self.really_all:
                self.examples = True
                self.modules = True
            if self.really_all:
                self.viewers = True

            # If we drop setuptools < 18.0, the remaining lines can probably be removed
            
            self.test_args = list(self._test_args())

            if noSuiteOrModule:
                # setuptools completely changed how it uses test_suite and test_args with v. 18.0
                # we do our best to keep it confused
                self.test_suite = "dummy"

        def _test_args(self):
            # can't seem to delegate a generator until Python 3.3
            
            if self.verbose:
                yield '--verbose'
            if self.test_suite:
                yield self.test_suite
                
            if self.viewers:
                print "*" * 60
                print "*" + "".center(58) + "*"
                print "*" + "ATTENTION".center(58) + "*"
                print "*" + "".center(58) + "*"
                print "*" + "Some of the following tests require user interaction".center(58) + "*"
                print "*" + "".center(58) + "*"
                print "*" * 60
                
                yield "fipy.viewers.testinteractive._suite"
            if self.modules:
                yield "fipy.testFiPy._suite"
            if self.examples:
                yield "examples.test._suite"

        def printPackageInfo(self):
            
            for pkg in ['fipy', 'numpy', 'pysparse', 'scipy', 'matplotlib', 'mpi4py']:

                try:
                    mod = __import__(pkg)

                    if hasattr(mod, '__version__'):
                        print pkg,'version',mod.__version__
                    else:
                        print pkg,'version not available'

                except ImportError, e:
                    print pkg,'is not installed'

                except Exception, e:
                    print pkg, 'version check failed:', e

            ## PyTrilinos
            try:
                import PyTrilinos
                print PyTrilinos.version()
            except ImportError, e:
                print pkg,'is not installed'
            except Exception, e:
                print pkg, 'version check failed:', e

            ## Mayavi uses a non-standard approach for storing its version nummber.
            try:
                from mayavi.__version__ import __version__ as mayaviversion
                print 'mayavi version', mayaviversion
            except ImportError, e:
                try:
                    from enthought.mayavi.__version__ import __version__ as mayaviversion
                    print 'enthought.mayavi version', mayaviversion
                except ImportError, e:
                    print 'enthought.mayavi is not installed'
                except Exception, e:
                    print 'enthought.mayavi version check failed:', e
            except Exception, e:
                print 'mayavi version check failed:', e

            ## Gmsh version
            try:
                from fipy.meshes.gmshMesh import gmshVersion
                gmshversion = gmshVersion()
                if gmshversion is None:
                    print 'gmsh is not installed'
                else:
                    print 'gmsh version',gmshversion
            except Exception, e:
                print 'gmsh version check failed:', e

        def run_tests(self):
            import sys
            if self.Trilinos or self.trilinos or self.no_pysparse:
                try:
                    ## The import scipy statement is added to allow
                    ## the --Trilinos tests to run without throwing a
                    ## segmentation fault. This is caused by weird
                    ## behavior in scipy and PyTrilinos depending on
                    ## the order in which modules are imported
                    try:
                        import scipy
                    except:
                        pass
                    import PyTrilinos
                except ImportError, a:
                    print >>sys.stderr, "!!! Trilinos library is not installed"
                    return

            if self.inline:
                try:
                    import weave
                except ImportError, a:
                    print >>sys.stderr, "!!! weave library is not installed"
                    return

            if self.pythoncompiled is not None:
                import os
                os.environ['PYTHONCOMPILED'] = self.pythoncompiled

            self.printPackageInfo()

            from pkg_resources import EntryPoint
            import unittest
            loader_ep = EntryPoint.parse("x="+self.test_loader)
            loader_class = loader_ep.load(require=False)

            from fipy.tools import numerix
            printoptions = numerix.get_printoptions()
            if "legacy" in printoptions:
                numerix.set_printoptions(legacy="1.13")

            try:
                unittest.main(
                    None, None, [unittest.__file__]+self.test_args,
                    testLoader = loader_class()
                    )
            except SystemExit, exitErr:
                # unittest.main(..., exit=...) not available until Python 2.7
                from fipy.tests.doctestPlus import report_skips
                report_skips()
                if self.timetests is not None:
                    pass
                else:
                    raise

            if "legacy" in printoptions:
                numerix.set_printoptions(legacy=printoptions["legacy"])

            if self.timetests is not None:
                from fipy.tests.doctestPlus import _DocTestTimes
                import numpy
                _DocTestTimes = numpy.rec.fromrecords(_DocTestTimes, formats='f8,S255', names='time,test')
                _DocTestTimes.sort(order=('time', 'test'))
                numpy.savetxt(self.timetests, _DocTestTimes[::-1], fmt="%8.4f\t%s")

            raise exitErr
    return _test
