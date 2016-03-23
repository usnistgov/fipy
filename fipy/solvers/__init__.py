from fipy.tools.parser import _parseSolver
from fipy.tools  import parallelComm as _parallelComm

from fipy.solvers.solver import *
__all__ = list(solver.__all__)

solver = _parseSolver()

def _envSolver(solver):
    import os
    if solver is None and 'FIPY_SOLVERS' in os.environ:
        solver = os.environ['FIPY_SOLVERS'].lower()
    return solver

solver = _envSolver(solver)

class SerialSolverError(Exception):
    def __init__(self, solver):
        super(SerialSolverError, self).__init__(solver + ' does not run in parallel')

if solver == "pysparse":
    if _parallelComm.Nproc > 1:
        raise SerialSolverError('pysparse')
    from fipy.solvers.pysparse import *
    __all__.extend(pysparse.__all__)
    from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
    _MeshMatrix =  _PysparseMeshMatrix

elif solver == "trilinos":
    from fipy.solvers.trilinos import *
    __all__.extend(trilinos.__all__)

    try:
        from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
        _MeshMatrix =  _PysparseMeshMatrix
    except ImportError:
        from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
        _MeshMatrix =  _TrilinosMeshMatrix

elif solver == "scipy":
    if _parallelComm.Nproc > 1:
        raise SerialSolverError('scipy')
    from fipy.solvers.scipy import *
    __all__.extend(scipy.__all__)
    from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
    _MeshMatrix = _ScipyMeshMatrix

elif solver == "pyamg":
    if _parallelComm.Nproc > 1:
        raise SerialSolverError('pyamg')
    from fipy.solvers.pyAMG import *
    __all__.extend(pyAMG.__all__)
    from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
    _MeshMatrix = _ScipyMeshMatrix

elif solver == "no-pysparse":
    from fipy.solvers.trilinos import *
    __all__.extend(trilinos.__all__)
    from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
    _MeshMatrix =  _TrilinosMeshMatrix

elif solver is None:
    # If no argument or environment variable, try importing them and seeing
    # what works

    exceptions = []

    try:
        if _parallelComm.Nproc > 1:
            raise SerialSolverError('pysparse')
        from fipy.solvers.pysparse import *
        __all__.extend(pysparse.__all__)
        solver = "pysparse"
        from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
        _MeshMatrix =  _PysparseMeshMatrix

    except (ImportError, SerialSolverError) as inst:
        exceptions.append(inst)

        try:
            from fipy.solvers.trilinos import *
            __all__.extend(trilinos.__all__)
            solver = "trilinos"
            try:
                from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix
                _MeshMatrix =  _PysparseMeshMatrix
            except ImportError:
                solver = "no-pysparse"
                from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrix
                _MeshMatrix =  _TrilinosMeshMatrix
        except ImportError as inst:
            exceptions.append(inst)

            try:
                if _parallelComm.Nproc > 1:
                    raise SerialSolverError('pyamg')
                from fipy.solvers.pyAMG import *
                __all__.extend(pyAMG.__all__)
                solver = "pyamg"
                from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
                _MeshMatrix = _ScipyMeshMatrix
            except (ImportError, SerialSolverError) as inst:
                exceptions.append(inst)

                try:
                    if _parallelComm.Nproc > 1:
                        raise SerialSolverError('scipy')
                    from fipy.solvers.scipy import *
                    __all__.extend(scipy.__all__)
                    solver = "scipy"
                    from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
                    _MeshMatrix = _ScipyMeshMatrix
                except (ImportError, SerialSolverError) as inst:
                    exceptions.append(inst)
                    import warnings
                    warnings.warn("Could not import any solver package. If you are using Trilinos, make sure you have all of the necessary Trilinos packages installed - Epetra, EpetraExt, AztecOO, Amesos, ML, and IFPACK.")
                    for inst in exceptions:
                        warnings.warn(inst.__class__.__name__ + ': ' + inst.message)


else:
    raise ImportError('Unknown solver package %s' % solver)


from fipy.tests.doctestPlus import register_skipper

register_skipper(flag='PYSPARSE_SOLVER',
                 test=lambda: solver == 'pysparse',
                 why="the PySparse solvers are not being used.",
                 skipWarning=True)
