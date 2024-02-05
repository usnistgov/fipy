from __future__ import print_function
from __future__ import unicode_literals
from builtins import str
from setuptools.command.test import test as _test
from future.utils import text_to_native_str
from future.utils import string_types
import unittest
import warnings
import sys

__all__ = ["DeprecationErroringTestProgram", "test"]
__all__ = [text_to_native_str(n) for n in __all__]

def _nativize_all(t):
    def _nativize(s):
        if isinstance(s, string_types):
            s = text_to_native_str(s)
        return s

    return tuple([_nativize(s) for s in t])

class DeprecationErroringTestProgram(unittest.TestProgram):
    """`TestProgram` that overrides inability of standard
    `TestProgram` to throw errors on `DeprecationWarning`
    """
    def runTests(self):
        self.warnings = None

        with warnings.catch_warnings():
            warnings.simplefilter(action="error", category=DeprecationWarning)

            # Don't raise noisy errors in
            # tvtk/tvtk_classes.zip/tvtk_classes/abstract_array.py or
            # tvtk/tvtk_classes.zip/tvtk_classes/algorithm.py
            warnings.filterwarnings(action="default", category=DeprecationWarning,
                                    message="invalid escape sequence.*")

            # Don't raise errors in skfmm.pfmm
            # due to deprecation of np.int in NumPy 1.20
            warnings.filterwarnings(action="default", category=DeprecationWarning,
                                    message="`np\.int` is a deprecated alias for the builtin `int`.*",
                                    module="skfmm\.pfmm")

            # Don't raise errors in scikits.umfpack
            # due to deprecation of pkg_resources.declare_namespace
            # raised in Python 3.7, but not newer
            # (combination of PyTrilinos and Gmsh(?) forces this old configuration)
            warnings.filterwarnings(action="default", category=DeprecationWarning,
                                    message="Deprecated call to `pkg_resources\.declare_namespace\('scikits'\)`.*")

            super(DeprecationErroringTestProgram, self).runTests()

class test(_test):
    description = str(_test.description) + ", for FiPy and its examples"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = _test.user_options + [
        ('inline', None, "run FiPy with inline compilation enabled"),
        ('pythoncompiled=', None, "directory in which to put weave's work product"),
        ('Trilinos', None, "run FiPy using Trilinos solvers"),
        ('Pysparse', None, "run FiPy using Pysparse solvers (default)"),
        ('trilinos', None, "run FiPy using Trilinos solvers"),
        ('pysparse', None, "run FiPy using Pysparse solvers (default)"),
        ('scipy', None, "run FiPy using SciPy solvers"),
        ('Scipy', None, "run FiPy using SciPy solvers"),
        ('petsc', None, "run FiPy using PETSc solvers"),
        ('no-pysparse', None, "run FiPy without using the Pysparse solvers"),
        ('pyamg', None, "run FiPy without using the PyAMG solvers"),
        ('pyamgx', None, "run FiPy using the pyamgx solvers"),
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
        ('deprecation-errors', None, "raise Exceptions for all DeprecationWarnings"),
       ]
    user_options = [_nativize_all(u) for u in user_options]

    def initialize_options(self):
        _test.initialize_options(self)

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
        self.pyamgx = False
        self.scipy = False
        self.petsc = False
        self.timetests = None
        self.skfmm = False
        self.lsmlib = False
        self.deprecation_errors = False

    def finalize_options(self):
        noSuiteOrModule = (self.test_suite is None
                           and self.test_module is None)

        _test.finalize_options(self)

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
            print("*" * 60)
            print("*" + "".center(58) + "*")
            print("*" + "ATTENTION".center(58) + "*")
            print("*" + "".center(58) + "*")
            print("*" + "Some of the following tests require user interaction".center(58) + "*")
            print("*" + "".center(58) + "*")
            print("*" * 60)

            yield "fipy.viewers.testinteractive._suite"
        if self.modules:
            yield "fipy.testFiPy._suite"
        if self.examples:
            yield "examples.test._suite"

    def printPackageInfo(self):
        from fipy.tools.logging.environment import package_info

        packages = package_info()

        print()
        package_width = max(len(word) for word in packages.keys())
        for package, version in packages.items():
            print("  ".join((package.ljust(package_width), version)))
        print()

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
            except ImportError as a:
                print("!!! Trilinos library is not installed", file=sys.stderr)
                return

        if self.pyamgx:
            try:
                ## Unregister the function pyamgx.finalize
                ## from atexit. This prevents
                ## pyamgx from printing an error message
                ## about memory leaks and a dump of leaked memory.
                ## The memory leaks happen because
                ## the tests do not use the pyamgx solvers
                ## "cleanly", i.e., they do not use the
                ## `with` statement.
                import pyamgx
                import atexit
                if hasattr(atexit, 'unregister'):
                    atexit.unregister(pyamgx.finalize)
                else:
                    atexit._exithandlers.remove(
                        (pyamgx.finalize, (), {}))
            except ImportError as e:
                print("!!! pyamgx package is not installed", file=sys.stederr)
                return

        if self.inline:
            try:
                import weave
            except ImportError as a:
                print("!!! weave library is not installed", file=sys.stderr)
                return

        if self.pythoncompiled is not None:
            import os
            os.environ['PYTHONCOMPILED'] = self.pythoncompiled

        self.printPackageInfo()

        from pkg_resources import EntryPoint
        loader_ep = EntryPoint.parse("x="+self.test_loader)
        loader_class = loader_ep.load(require=False)

        from fipy.tools import numerix
        printoptions = numerix.get_printoptions()
        if "legacy" in printoptions:
            numerix.set_printoptions(legacy="1.13")

        if self.deprecation_errors:
            test_program_class = DeprecationErroringTestProgram
        else:
            test_program_class = unittest.TestProgram

        try:
            test_program_class(
                None, None, [unittest.__file__]+self.test_args,
                testLoader = loader_class()
                )
        except SystemExit as exitErr:
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

            # Only write on proc 0.
            # Doesn't use FiPy comms because this command can
            # be run outside of FiPy (`python setup.py test`)
            try:
                from mpi4py import MPI
                procID = MPI.COMM_WORLD.rank
                barrier = MPI.COMM_WORLD.barrier
            except:
                procID = 0
                def barrier(*args):
                    pass

            if procID == 0:
                numpy.savetxt(self.timetests, _DocTestTimes[::-1],
                              delimiter='\t',
                              header="time\tmodule", comments='',
                              fmt=("%.18e", "%s"))

            barrier()

        raise exitErr
