#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "trilinosNonlinearSolver.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: David Saylor <David.Saylor@fda.hhs.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

from PyTrilinos import Epetra
from PyTrilinos import NOX

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.tools import parallelComm

class _NOXInterface(NOX.Epetra.Interface.Required, NOX.Epetra.Interface.Jacobian):
    def __init__(self, solver):
        NOX.Epetra.Interface.Required.__init__(self)
        NOX.Epetra.Interface.Jacobian.__init__(self)
        self.solver = solver

    def solve(self, dt=None):
        self.dt = dt

        self.solver.equation._prepareLinearSystem(var=None, solver=self.solver, boundaryConditions=(), dt=1.)

        globalMatrix, nonOverlappingVector, nonOverlappingRHSvector, overlappingVector = self.solver._globalMatrixAndVectors

        self.colMap = globalMatrix.colMap
        self.domainMap = globalMatrix.domainMap

        if self.solver.jacobian is None:
            # Define the Jacobian interface/operator
            jacInterface = NOX.Epetra.MatrixFree(self.solver.nlParams["Printing"], self, nonOverlappingVector)
            jacobian = jacInterface
        else:
            jacInterface = self
            self.solver.jacobian.cacheMatrix()
            self.jacsolver = self.solver.jacobian._prepareLinearSystem(var=None, solver=_DummyJacobianSolver(), boundaryConditions=(), dt=1.)
            jacobian = (self.solver.jacobian.matrix.asTrilinosMeshMatrix()).matrix

        # Define the Preconditioner interface/operator
        fdc = NOX.Epetra.FiniteDifferenceColoring(self.solver.nlParams["Printing"], self,
                                                  nonOverlappingVector, globalMatrix.matrix.Graph(), True)

        noxSolver = NOX.Epetra.defaultSolver(initGuess=nonOverlappingVector,
                                             reqInterface=self,
                                             jacInterface=jacInterface, jacobian=jacobian,
                                             precInterface=fdc, preconditioner=fdc,
                                             nlParams=self.solver.nlParams,
                                             absTol=self.solver.tolerance, relTol=0.5, maxIters=self.solver.iterations,
                                             updateTol=None, wAbsTol=None, wRelTol=None)

        status = noxSolver.solve()

        self.solver._deleteGlobalMatrixAndVectors()

        return status


    def computeJacobian(self, u, Jac):
        try:

            overlappingVector = Epetra.Vector(self.colMap, self.solver.var)

            overlappingVector.Import(u,
                                     Epetra.Import(self.colMap,
                                                   self.domainMap),
                                     Epetra.Insert)

            self.solver.var.value = overlappingVector

            self.solver.jacobian.matrix.trilinosMatrix = Jac
            self.solver.jacobian._prepareLinearSystem(var=None, solver=self.jacsolver, boundaryConditions=(), dt=1.)

            return True

        except Exception, e:
            print "TrilinosNonlinearSolver.computeJacobian() has thrown an exception:"
            print str(type(e))[18:-2] + ":", e
            return False


    def computeF(self, u, F, flag):
        try:
            overlappingVector = Epetra.Vector(self.colMap, self.solver.var)

            overlappingVector.Import(u,
                                     Epetra.Import(self.colMap,
                                                   self.domainMap),
                                     Epetra.Insert)

            self.solver.var.value = overlappingVector
            F[:] = self.solver.equation.justResidualVector(dt=self.dt)

            return True

        except Exception, e:
            print "TrilinosNonlinearSolver.computeF() has thrown an exception:"
            print str(type(e))[18:-2] + ":", e
            return False

class _DummyJacobianSolver(TrilinosSolver):
    pass

__all__ = ["TrilinosNonlinearSolver"]

class TrilinosNonlinearSolver(TrilinosSolver):
    def __init__(self, equation, jacobian=None, tolerance=1e-10, iterations=1000,
                 printingOptions=None, solverOptions=None, linearSolverOptions=None,
                 lineSearchOptions=None, directionOptions=None, newtonOptions=None):
        TrilinosSolver.__init__(self, tolerance=tolerance, iterations=iterations, precon=None)

        self.equation = equation
        self.jacobian = jacobian

        self.nlParams = NOX.Epetra.defaultNonlinearParameters(parallelComm.epetra_comm, 2)
        self.nlParams["Printing"] = printingOptions or {
            'Output Precision': 3,
            'MyPID': 0,
            'Output Information': NOX.Utils.OuterIteration,
            'Output Processor': 0
        }
        self.nlParams["Solver Options"] =  solverOptions or {
            'Status Test Check Type': 'Complete',
            'Rescue Bad Newton Solve': 'True'
        }
        self.nlParams["Linear Solver"] = linearSolverOptions or {
            'Aztec Solver': 'GMRES',
            'Tolerance': 0.0001,
            'Max Age Of Prec': 5,
            'Max Iterations': 20,
            'Preconditioner': 'Ifpack'
        }
        self.nlParams["Line Search"] = lineSearchOptions or {'Method': "Polynomial"}
        self.nlParams["Direction"] = directionOptions or {'Method': 'Newton'}
        self.nlParams["Newton"] = newtonOptions or {'Forcing Term Method': 'Type 2'}

        self.nox = _NOXInterface(solver=self)

    def solve(self, dt=None):
        output = self.nox.solve(dt=dt)

#         if 'FIPY_VERBOSE_SOLVER' in os.environ:
#             status = Solver.GetAztecStatus()
#
#             from fipy.tools.debug import PRINT
#             PRINT('iterations: %d / %d' % (status[AztecOO.AZ_its], self.iterations))
#             failure = {AztecOO.AZ_normal : 'AztecOO.AZ_normal',
#                        AztecOO.AZ_param : 'AztecOO.AZ_param',
#                        AztecOO.AZ_breakdown : 'AztecOO.AZ_breakdown',
#                        AztecOO.AZ_loss : 'AztecOO.AZ_loss',
#                        AztecOO.AZ_ill_cond : 'AztecOO.AZ_ill_cond',
#                        AztecOO.AZ_maxits : 'AztecOO.AZ_maxits'}
#
#             PRINT('failure',failure[status[AztecOO.AZ_why]])
#
#             PRINT('AztecOO.AZ_r:',status[AztecOO.AZ_r])
#             PRINT('AztecOO.AZ_scaled_r:',status[AztecOO.AZ_scaled_r])
#             PRINT('AztecOO.AZ_rec_r:',status[AztecOO.AZ_rec_r])
#             PRINT('AztecOO.AZ_solve_time:',status[AztecOO.AZ_solve_time])
#             PRINT('AztecOO.AZ_Aztec_version:',status[AztecOO.AZ_Aztec_version])

        return output
