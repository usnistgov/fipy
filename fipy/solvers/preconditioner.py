from builtins import object
__docformat__ = 'restructuredtext'

__all__ = ["Preconditioner",
           "SolverModifyingPreconditioner",
           "MatrixModifyingPreconditioner"]

class Preconditioner(object):
    """Base class for solver preconditioners.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    pass


class SolverModifyingPreconditioner(Preconditioner):
    """Base class for preconditioners that modify a :class:`~fipy.solvers.solver.Solver`.
    """

    def _applyToSolver(self, solver, matrix):
        """Modify `solver` to apply preconditioning to `matrix`.

        Parameters
        ----------
        solver
            The solver to modify with preconditioner.
        matrix
            The matrix the preconditioner applies to.

        Returns
        -------
        None
        """
        raise NotImplementedError


class MatrixModifyingPreconditioner(Preconditioner):
    """Base class for preconditioners that modify a :class:`~fipy.matrices.sparseMatrix._SparseMatrix`.
    """

    def _applyToMatrix(self, matrix):
        """Create a preconditioner for `matrix`.

        Returns
        -------
        preconditioner : object
            Preconditioning object appropriate for this solver suite.
        matrix : :class:`~fipy.matrices.sparseMatrix._SparseMatrix`
            Matrix, possibly restructured to facilitate applying
            preconditioner.
        """
        raise NotImplementedError

