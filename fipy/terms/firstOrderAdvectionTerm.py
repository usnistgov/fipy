from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = ['FirstOrderAdvectionTerm']
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from fipy.tools import numerix
from fipy.tools.numerix import MA

from fipy.terms.nonDiffusionTerm import _NonDiffusionTerm

class FirstOrderAdvectionTerm(_NonDiffusionTerm):
    r"""

    The `FirstOrderAdvectionTerm` object constructs the b vector contribution
    for the advection term given by

    .. math::

       u \abs{\nabla \phi}

    from the advection equation given by:

    .. math::

       \frac{\partial \phi}{\partial t} + u \abs{\nabla \phi} = 0

    The construction of the gradient magnitude term requires upwinding.
    The formula used here is given by:

    .. math::

       u_P \abs{\nabla \phi}_P = \max \left( u_P , 0 \right) \left[  \sum_A \min \left( \frac{ \phi_A - \phi_P } { d_{AP}}, 0 \right)^2 \right]^{1/2} +  \min \left( u_P , 0 \right) \left[  \sum_A \max \left( \frac{ \phi_A - \phi_P } { d_{AP}}, 0 \right)^2 \right]^{1/2}

    Here are some simple test cases for this problem:

    >>> from fipy.meshes import Grid1D
    >>> from fipy.solvers import *
    >>> SparseMatrix = LinearLUSolver()._matrixClass
    >>> mesh = Grid1D(dx = 1., nx = 3)
    >>> from fipy.variables.cellVariable import CellVariable

    Trivial test:

    >>> var = CellVariable(value = numerix.zeros(3, 'd'), mesh = mesh)
    >>> v, L, b = FirstOrderAdvectionTerm(0.)._buildMatrix(var, SparseMatrix)
    >>> print(numerix.allclose(b, numerix.zeros(3, 'd'), atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    Less trivial test:

    >>> var = CellVariable(value = numerix.arange(3), mesh = mesh)
    >>> v, L, b = FirstOrderAdvectionTerm(1.)._buildMatrix(var, SparseMatrix)
    >>> print(numerix.allclose(b, numerix.array((0., -1., -1.)), atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    Even less trivial

    >>> var = CellVariable(value = numerix.arange(3), mesh = mesh)
    >>> v, L, b = FirstOrderAdvectionTerm(-1.)._buildMatrix(var, SparseMatrix)
    >>> print(numerix.allclose(b, numerix.array((1., 1., 0.)), atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    Another trivial test case (more trivial than a trivial test case
    standing on a harpsichord singing "trivial test cases are here again")

    >>> vel = numerix.array((-1, 2, -3))
    >>> var = CellVariable(value = numerix.array((4, 6, 1)), mesh = mesh)
    >>> v, L, b = FirstOrderAdvectionTerm(vel)._buildMatrix(var, SparseMatrix)
    >>> print(numerix.allclose(b, -vel * numerix.array((2, numerix.sqrt(5**2 + 2**2), 5)), atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    Somewhat less trivial test case:

    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
    >>> vel = numerix.array((3, -5, -6, -3))
    >>> var = CellVariable(value = numerix.array((3, 1, 6, 7)), mesh = mesh)
    >>> v, L, b = FirstOrderAdvectionTerm(vel)._buildMatrix(var, SparseMatrix)
    >>> answer = -vel * numerix.array((2, numerix.sqrt(2**2 + 6**2), 1, 0))
    >>> print(numerix.allclose(b, answer, atol = 1e-10)) # doctest: +PROCESSOR_0
    True
    """

    def __init__(self, coeff = None):
        _NonDiffusionTerm.__init__(self)
        self.geomCoeff = coeff

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, equation=None, transientGeomCoeff=None, diffusionGeomCoeff=None):

        oldArray = var.old

        mesh = var.mesh
        NCells = mesh.numberOfCells
        NCellFaces = mesh._maxFacesPerCell

        cellValues = numerix.repeat(oldArray[numerix.newaxis, ...], NCellFaces, axis = 0)

        cellIDs = numerix.repeat(numerix.arange(NCells)[numerix.newaxis, ...], NCellFaces, axis = 0)
        cellToCellIDs = mesh._cellToCellIDs

        if NCells > 0:
            cellToCellIDs = MA.where(MA.getmask(cellToCellIDs), cellIDs, cellToCellIDs)

            adjacentValues = numerix.take(oldArray, cellToCellIDs)

            differences = self._getDifferences(adjacentValues, cellValues, oldArray, cellToCellIDs, mesh)
            differences = MA.filled(differences, 0)

            minsq = numerix.sqrt(numerix.sum(numerix.minimum(differences, numerix.zeros((NCellFaces, NCells), 'l'))**2, axis=0))
            maxsq = numerix.sqrt(numerix.sum(numerix.maximum(differences, numerix.zeros((NCellFaces, NCells), 'l'))**2, axis=0))

            coeff = numerix.array(self._getGeomCoeff(var))

            coeffXdifferences = coeff * ((coeff > 0.) * minsq + (coeff < 0.) * maxsq)
        else:
            coeffXdifferences = 0.

        return (var, SparseMatrix(mesh=var.mesh), -coeffXdifferences * mesh.cellVolumes)

    def _getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs, mesh):
        return (adjacentValues - cellValues) / mesh._cellToCellDistances

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        solver = solver or super(FirstOrderAdvectionTerm, self)._getDefaultSolver(var, solver, *args, **kwargs)

        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve asymmetric matrices" % solver)

        import fipy.solvers.solver
        if fipy.solvers.solver_suite in ['trilinos', 'no-pysparse']:
            from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner
            from fipy.solvers.trilinos.linearGMRESSolver import LinearGMRESSolver
            return solver or LinearGMRESSolver(precon=JacobiPreconditioner(), *args, **kwargs)
        elif fipy.solvers.solver_suite == 'pyamg':
            from fipy.solvers.pyAMG.linearGeneralSolver import LinearGeneralSolver
            return solver or LinearGeneralSolver(tolerance=1e-15, iterations=2000, *args, **kwargs)
        else:
            from fipy.solvers import DefaultAsymmetricSolver
            return solver or DefaultAsymmetricSolver(*args, **kwargs)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


