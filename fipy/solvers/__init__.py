from __future__ import unicode_literals
from builtins import str

import logging

_log = logging.getLogger(__name__)

import os
from importlib import import_module

from fipy.tools.parser import _parseSolver

from fipy.solvers.solver import *

_desired_solver = _parseSolver()

if _desired_solver is None and 'FIPY_SOLVERS' in os.environ:
    _desired_solver = os.environ['FIPY_SOLVERS'].lower()
del os
    
try:
    from mpi4py import MPI
    _Nproc = MPI.COMM_WORLD.size
    del MPI
except ImportError:
    _Nproc = 1

_exceptions = []

class SerialSolverError(Exception):
    def __init__(self):
        super(SerialSolverError, self).__init__('solver does not run in parallel')

_exceptions = {}

def _import_matrix(suite, kind):
    """`from fipy.matrices.suiteMatrix import _SuiteKindMatrix`
    """
    m = import_module("fipy.matrices.{}Matrix".format(suite.lower()))
    return getattr(m, "_{}{}Matrix".format(suite, kind))

def _import_mesh_matrices(suite):
    """Import row, column, and general mesh matrices from `suite`
    """
    _RowMeshMatrix = _import_matrix(suite=suite, kind="RowMesh")
    _ColMeshMatrix = _import_matrix(suite=suite, kind="ColMesh")
    _MeshMatrix = _import_matrix(suite=suite, kind="Mesh")

    return _RowMeshMatrix, _ColMeshMatrix, _MeshMatrix

solver = None

from fipy.tools.comms.dummyComm import DummyComm
serialComm, parallelComm = DummyComm(), DummyComm()

if solver is None and _desired_solver in ["pysparse", None]:
    try:
        if _Nproc > 1:
            raise SerialSolverError()
        from fipy.solvers.pysparse import *
        _mesh_matrices = _import_mesh_matrices(suite="Pysparse")
        solver = "pysparse"
    except Exception as inst:
        _exceptions["pysparse"] = inst

if solver is None and _desired_solver in ["petsc", None]:
    try:
        import petsc4py
        petsc4py.init()

        from fipy.solvers.petsc import *

        from fipy.solvers.petsc.comms.serialPETScCommWrapper import SerialPETScCommWrapper
        serialComm = SerialPETScCommWrapper()

        if _Nproc > 1:
            from fipy.solvers.petsc.comms.parallelPETScCommWrapper import ParallelPETScCommWrapper
            parallelComm = ParallelPETScCommWrapper()
        else:
            parallelComm = SerialPETScCommWrapper()

        _mesh_matrices = _import_mesh_matrices(suite="PETSc")
        solver = "petsc"
    except Exception as inst:
        _exceptions["petsc"] = inst

if solver is None and _desired_solver in ["trilinos", "no-pysparse", None]:
    try:
        from fipy.solvers.trilinos import *
        
        from fipy.solvers.trilinos.comms.serialEpetraCommWrapper import SerialEpetraCommWrapper
        serialComm = SerialEpetraCommWrapper()

        if _Nproc > 1:
            from fipy.solvers.trilinos.comms.parallelEpetraCommWrapper import ParallelEpetraCommWrapper
            parallelComm = ParallelEpetraCommWrapper()
        else:
            parallelComm = SerialEpetraCommWrapper()

        if _desired_solver != "no-pysparse":
            try:
                _mesh_matrices = _import_mesh_matrices(suite="Pysparse")
                solver = "trilinos"
            except ImportError:
                pass
                
        if solver is None:
            # no-pysparse requested or pysparseMatrix failed to import
            _mesh_matrices = _import_mesh_matrices(suite="Trilinos")
            solver = "no-pysparse"
    except Exception as inst:
        _exceptions["trilinos"] = inst

if solver is None and _desired_solver in ["scipy", None]:
    try:
        if _Nproc > 1:
            raise SerialSolverError()
        from fipy.solvers.scipy import *
        _mesh_matrices = _import_mesh_matrices(suite="Scipy")
        solver = "scipy"
    except Exception as inst:
        _exceptions["scipy"] = inst

if solver is None and _desired_solver in ["pyamg", None]:
    try:
        if _Nproc > 1:
            raise SerialSolverError()
        from fipy.solvers.pyAMG import *
        _mesh_matrices = _import_mesh_matrices(suite="Scipy")
        solver = "pyamg"
    except Exception as inst:
        _exceptions["pyamg"] = inst

if solver is None and _desired_solver in ["pyamgx", None]:
    try:
        if _parallelComm.Nproc > 1:
            raise  SerialSolverError('pyamgx')
        from fipy.solvers.pyamgx import *
        _mesh_matrices = _import_mesh_matrices(suite="Scipy")
        solver = "pyamgx"
    except Exception as inst:
        _exceptions["pyamgx"] = inst

if solver is None:
    if _desired_solver is None:
        raise ImportError('Unable to load a solver: %s' % str(_exceptions))
    else:
        if len(_exceptions) > 0:
            raise ImportError('Unable to load solver %s: %s' % (_desired_solver, _exceptions[_desired_solver]))
        else:
            raise ImportError('Unknown solver package %s' % _desired_solver)

# don't unpack until here in order to keep code above more succinct
_RowMeshMatrix, _ColMeshMatrix, _MeshMatrix = _mesh_matrices

from fipy.tests.doctestPlus import register_skipper

register_skipper(flag='PYSPARSE_SOLVER',
                 test=lambda: solver == 'pysparse',
                 why="the Pysparse solvers are not being used.",
                 skipWarning=True)

register_skipper(flag='NOT_PYAMGX_SOLVER',
                 test=lambda: solver != 'pyamgx',
                 why="the PyAMGX solver is being used.",
                 skipWarning=True)
del register_skipper

_log.info("Solver suite is %s", solver)
_log.debug("DefaultSolver is %s", DefaultSolver)
_log.debug("DefaultAsymmetricSolver is %s", DefaultAsymmetricSolver)
_log.debug("DummySolver is %s", DummySolver)
_log.debug("GeneralSolver is %s", GeneralSolver)
_log.debug("serialComm is %s", serialComm)
_log.debug("parallelComm is %s", parallelComm)
