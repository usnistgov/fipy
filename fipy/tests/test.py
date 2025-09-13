import contextlib
import os
from pathlib import Path
import sys
import typer
from typer import Argument, Option
from typing_extensions import Annotated
import unittest
import warnings

__all__ = ["app"]

## cribbed from obsolete setuptools.command.test

@contextlib.contextmanager
def project_on_sys_path(include_dists=[]):
    old_path = sys.path[:]
    old_modules = sys.modules.copy()

    try:
        project_path = Path(".").absolute().as_posix()
        sys.path.insert(0, project_path)
        with paths_on_pythonpath([project_path]):
            yield
    finally:
        sys.path[:] = old_path
        sys.modules.clear()
        sys.modules.update(old_modules)

@contextlib.contextmanager
def paths_on_pythonpath(paths):
    """
    Add the indicated paths to the head of the :envvar:`PYTHONPATH`
    environment variable so that subprocesses will also see the packages at
    these paths.

    Do this in a context that restores the value on exit.
    """
    nothing = object()
    orig_pythonpath = os.environ.get('PYTHONPATH', nothing)
    current_pythonpath = os.environ.get('PYTHONPATH', '')
    try:
        prefix = os.pathsep.join(paths)
        to_join = filter(None, [prefix, current_pythonpath])
        new_path = os.pathsep.join(to_join)
        if new_path:
            os.environ['PYTHONPATH'] = new_path

        yield
    finally:
        if orig_pythonpath is nothing:
            os.environ.pop('PYTHONPATH', None)
        else:
            os.environ['PYTHONPATH'] = orig_pythonpath

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
                                    message=r"`np\.int` is a deprecated alias for the builtin `int`.*",
                                    module=r"skfmm\.pfmm")

            # Don't raise errors in scikits.umfpack
            # due to deprecation of pkg_resources.declare_namespace
            # raised in Python 3.7, but not newer
            # (combination of PyTrilinos and Gmsh(?) forces this old configuration)
            warnings.filterwarnings(action="default", category=DeprecationWarning,
                                    message=r"Deprecated call to `pkg_resources\.declare_namespace\('scikits'\)`.*")

            super(DeprecationErroringTestProgram, self).runTests()

def _test_args(verbose, module_or_suite, viewers, modules, examples):
    # can't seem to delegate a generator until Python 3.3

    if verbose:
        yield '--verbose'
    if module_or_suite:
        yield module_or_suite

    if viewers:
        print("*" * 60)
        print("*" + "".center(58) + "*")
        print("*" + "ATTENTION".center(58) + "*")
        print("*" + "".center(58) + "*")
        print("*" + "Some of the following tests require user interaction".center(58) + "*")
        print("*" + "".center(58) + "*")
        print("*" * 60)

        yield "fipy.viewers.testinteractive._suite"
    if modules:
        yield "fipy.testFiPy._suite"
    if examples:
        yield "examples.test._suite"

def _printPackageInfo():
    from fipy.tools.logging.environment import package_info

    packages = package_info()

    print()
    package_width = max(len(word) for word in packages.keys())
    for package, version in packages.items():
        print("  ".join((package.ljust(package_width), version)))
    print()

def _initialize_trilinos():
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

def _initialize_pyamgx():
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

def _initialize_weave():
    try:
        import weave
    except ImportError as a:
        print("!!! weave library is not installed", file=sys.stderr)
        return

def run_tests(test_program_class, test_args, timetests=None):
    _printPackageInfo()

    from fipy.tools import numerix
    printoptions = numerix.get_printoptions()
    if "legacy" in printoptions:
        numerix.set_printoptions(legacy="1.13")

    try:
        test_program_class(
            None, None, [unittest.__file__]+test_args,
            testLoader = unittest.TestLoader()
            )
    except SystemExit as exitErr:
        # unittest.main(..., exit=...) not available until Python 2.7
        from fipy.tests.doctestPlus import report_skips
        report_skips()
        if timetests is not None:
            pass
        else:
            raise

    if "legacy" in printoptions:
        numerix.set_printoptions(legacy=printoptions["legacy"])

    if timetests is not None:
        from fipy.tests.doctestPlus import _DocTestTimes
        import numpy
        _DocTestTimes = numpy.rec.fromrecords(_DocTestTimes, formats='f8,S255', names='time,test')
        _DocTestTimes.sort(order=('time', 'test'))

        # Only write on proc 0.
        # Doesn't use FiPy comms because this command can
        # be run outside of FiPy (`fipy_test`)
        try:
            from mpi4py import MPI
            procID = MPI.COMM_WORLD.rank
            barrier = MPI.COMM_WORLD.barrier
        except:
            procID = 0
            def barrier(*args):
                pass

        if procID == 0:
            numpy.savetxt(timetests, _DocTestTimes[::-1],
                          delimiter='\t',
                          header="time\tmodule", comments='',
                          fmt=("%.18e", "%s"))

        barrier()

        if 'exitErr' in locals():
            raise exitErr

def main(
    module_or_suite: Annotated[
        str,
        Argument(help="module or test suite to test")
    ] = None,
    inline: Annotated[
        bool,
        Option(help="run FiPy with inline compilation enabled")
    ] = False,
    pythoncompiled: Annotated[
        Path,
        Option(help="directory in which to put weave's work product")
    ] = None,
    trilinos: Annotated[
        bool,
        Option(help="run FiPy using Trilinos solvers")
    ] = False,
    scipy: Annotated[
        bool,
        Option(help="run FiPy using SciPy solvers")
    ] = False,
    petsc: Annotated[
        bool,
        Option(help="run FiPy using PETSc solvers")
    ] = False,
    pyamgx: Annotated[
        bool,
        Option(help="run FiPy using pyamgx solvers")
    ] = False,
    all: Annotated[
        bool,
        Option(help="run all non-interactive FiPy tests (default)")
    ] = None,
    really_all: Annotated[
        bool,
        Option(help="run *all* FiPy tests (including those requiring user input)")
    ] = False,
    examples: Annotated[
        bool,
        Option(help="test FiPy examples")
    ] = False,
    modules: Annotated[
        bool,
        Option(help="test FiPy code modules")
    ] = False,
    viewers: Annotated[
        bool,
        Option(help="test FiPy viewer modules (requires user input)")
    ] = False,
    cache: Annotated[
        bool,
        Option(help="run FiPy with Variable caching")
    ] = False,
    no_cache: Annotated[
        bool,
        Option(help="run FiPy without Variable caching")
    ] = True,
    timetests: Annotated[
        Path,
        Option(help="file in which to put time spent on each test")
    ] = None,
    skfmm: Annotated[
        bool,
        Option(help="run FiPy using the Scikit-fmm level set solver")
    ] = False,
    lsmlib: Annotated[
        bool,
        Option(help="run FiPy using the LSMLIB level set solver")
    ] = False,
    deprecation_errors: Annotated[
        bool,
        Option(help="raise Exceptions for all DeprecationWarnings")
    ] = False,
    verbose: Annotated[
        bool,
        Option(help="blah blah verbose")
    ] = True
):
    if trilinos:
        _initialize_trilinos()

    if pyamgx:
        _initialize_pyamgx()

    if inline:
        _initialize_weave()

    if pythoncompiled is not None:
        os.environ['PYTHONCOMPILED'] = pythoncompiled.absolute().as_posix()

    if not (examples or modules or viewers or module_or_suite):
        all = True
    if all or really_all:
        examples = True
        modules = True
    if really_all:
        viewers = True

    test_args = list(_test_args(verbose, module_or_suite, viewers, modules, examples))

    if deprecation_errors:
        test_program_class = DeprecationErroringTestProgram
    else:
        test_program_class = unittest.TestProgram

    # Run tests with current working directory on path
    # so that examples can be found
    with project_on_sys_path():
        run_tests(test_program_class=test_program_class,
                  test_args=test_args,
                  timetests=timetests)

app = typer.Typer()
app.command()(main)

if __name__ == "__main__":
    app()
