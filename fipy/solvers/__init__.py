from __future__ import unicode_literals
from fipy.tools.parser import _parseSolver
from fipy.tools  import parallelComm as _parallelComm

from fipy.solvers.solver import *
__all__ = list(solver.__all__)
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

_desired_solver = _parseSolver()

import os
if _desired_solver is None and 'FIPY_SOLVERS' in os.environ:
    _desired_solver = os.environ['FIPY_SOLVERS'].lower()
    
exceptions = []

class SerialSolverError(Exception):
    def __init__(self, solver):
        super(SerialSolverError, self).__init__(solver + ' does not run in parallel')

solver = None

if solver is None and _desired_solver in ["pysparse", None]:
    try:
        if _parallelComm.Nproc > 1:
            raise SerialSolverError('pysparse')
        from fipy.solvers.pysparse import *
        __all__.extend(pysparse.__all__)
        from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
        _MeshMatrix =  _PysparseMeshMatrix
        solver = "pysparse"
    except Exception as inst:
        exceptions.append(inst)

if solver is None and _desired_solver in ["petsc", None]:
    try:
        from fipy.solvers.petsc import *
        __all__.extend(petsc.__all__)

        from fipy.matrices.petscMatrix import _PETScMeshMatrix
        _MeshMatrix =  _PETScMeshMatrix
        solver = "petsc"
    except Exception as inst:
        exceptions.append(inst)

if solver is None and _desired_solver in ["trilinos", "no-pysparse", None]:
    try:
        from fipy.solvers.trilinos import *
        __all__.extend(trilinos.__all__)

        if _desired_solver != "no-pysparse":
            try:
                from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
                _MeshMatrix =  _PysparseMeshMatrix
                solver = "trilinos"
            except ImportError:
                pass
                
        if solver is None:
            # no-pysparse requested or pysparseMatrix failed to import
            from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
            _MeshMatrix =  _TrilinosMeshMatrix
            solver = "no-pysparse"
    except Exception as inst:
        exceptions.append(inst)

if solver is None and _desired_solver in ["scipy", None]:
    try:
        if _parallelComm.Nproc > 1:
            raise  SerialSolverError('scipy')
        from fipy.solvers.scipy import *
        __all__.extend(scipy.__all__)
        from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
        _MeshMatrix = _ScipyMeshMatrix
        solver = "scipy"
    except Exception as inst:
        exceptions.append(inst)

if solver is None and _desired_solver in ["pyamg", None]:
    try:
        if _parallelComm.Nproc > 1:
            raise  SerialSolverError('pyamg')
        from fipy.solvers.pyAMG import *
        __all__.extend(pyAMG.__all__)
        from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
        _MeshMatrix = _ScipyMeshMatrix
        solver = "pyamg"
    except Exception as inst:
        exceptions.append(inst)

if solver is None and _desired_solver in ["pyamgx", None]:
    try:
        if _parallelComm.Nproc > 1:
            raise  SerialSolverError('pyamgx')
        from fipy.solvers.pyamgx import *
        __all__.extend(pyamgx.__all__)
        from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
        _MeshMatrix = _ScipyMeshMatrix
        solver = "pyamgx"
    except Exception as inst:
        exceptions.append(inst)

if solver is None:
    if _desired_solver is None:
        raise ImportError('Unable to load a solver: %s' % [str(e) for e in exceptions])
    else:
        if len(exceptions) > 0:
            raise ImportError('Unable to load solver %s: %s' % (_desired_solver, [str(e) for e in exceptions]))
        else:
            raise ImportError('Unknown solver package %s' % _desired_solver)

from fipy.tests.doctestPlus import register_skipper

register_skipper(flag='PYSPARSE_SOLVER',
                 test=lambda: solver == 'pysparse',
                 why="the Pysparse solvers are not being used.",
                 skipWarning=True)

register_skipper(flag='NOT_PYAMGX_SOLVER',
                 test=lambda: solver != 'pyamgx',
                 why="the PyAMGX solver is being used.",
                 skipWarning=True)
