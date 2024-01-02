from __future__ import unicode_literals

from scipy.sparse import csr_matrix, linalg

import pyamgx

from fipy.solvers.solver import Solver
from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.tools import numerix

__all__ = ["PyAMGXSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PyAMGXSolver(Solver):

    # AMGX configuration options
    CONFIG_DICT = {}

    #: Default smoother to apply to the ???
    DEFAULT_SMOOTHER = None

    def __init__(self, tolerance="default", criterion="default",
                 iterations="default", precon="default", smoother="default", **kwargs):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pyamgx.preconditioners.PyAMGXPreconditioner, optional
        smoother : ~fipy.solvers.pyamgx.smoothers.Smoother, optional
        **kwargs
            Other AMGX solver options
        """
        super(PyAMGXSolver, self).__init__(tolerance=tolerance, criterion=criterion, iterations=iterations)

        # update solver config:
        self.config_dict = self.CONFIG_DICT.copy()

        self.config_dict["solver"]["max_iters"] = self.iterations

        if self.precon is not None:
            self.precon._applyToSolver(self.config_dict["solver"])

        smoother = self.value_or_default(smoother, self.default_smoother)
        if smoother is not None:
            smoother._applyToSolver(self.config_dict["solver"])

        self.config_dict["solver"].update(kwargs)

        # create AMGX objects:
        self.cfg = pyamgx.Config().create_from_dict(self.config_dict)
        self.resources = pyamgx.Resources().create_simple(self.cfg)
        self.x_gpu = pyamgx.Vector().create(self.resources)
        self.b_gpu = pyamgx.Vector().create(self.resources)
        self.A_gpu = pyamgx.Matrix().create(self.resources)

    @property
    def default_smoother(self):
        if self.DEFAULT_SMOOTHER is not None:
            # instantiate DEFAULT_SMOOTHER class
            return self.DEFAULT_SMOOTHER()
        else:
            return None

    def _destroy_AMGX(self):
        # destroy AMGX objects:
        # self.resources apparently doesn't need to be destroyed
        if hasattr(self, "A_gpu"):
            self.A_gpu.destroy()
            del self.A_gpu
        if hasattr(self, "b_gpu"):
            self.b_gpu.destroy()
            del self.b_gpu
        if hasattr(self, "x_gpu"):
            self.x_gpu.destroy()
            del self.x_gpu
        if hasattr(self, "cfg"):
            self.cfg.destroy()
            del self.cfg

    def __exit__(self, *args):
        self._destroy_AMGX()

    def __del__(self):
        self._destroy_AMGX()

    @property
    def _matrixClass(self):
        return _ScipyMeshMatrix

    def _rhsNorm(self, L, x, b):
        return numerix.L2norm(b)

    def _matrixNorm(self, L, x, b):
        return linalg.norm(L.matrix, ord=numerix.inf)

    def _adaptLegacyTolerance(self, L, x, b):
        return self._adaptInitialTolerance(L, x, b)

    def _adaptUnscaledTolerance(self, L, x, b):
        return (1., "ABSOLUTE")

    def _adaptRHSTolerance(self, L, x, b):
        return (self._rhsNorm(L, x, b), "ABSOLUTE")

    def _adaptMatrixTolerance(self, L, x, b):
        return (self._matrixNorm(L, x, b), "ABSOLUTE")

    def _adaptInitialTolerance(self, L, x, b):
        return (1., "RELATIVE_INI_CORE")

    def _solve_(self, L, x, b):
        # transfer data from CPU to GPU
        self.x_gpu.upload(x)
        self.b_gpu.upload(b)
        self.A_gpu.upload_CSR(L.matrix)

        tolerance_scale, suite_criterion = self._adaptTolerance(L, x, b)
        config_dict = self.config_dict.copy()
        config_dict["solver"]["monitor_residual"] = 1
        config_dict["solver"]["store_res_history"] = 1
        config_dict["solver"]["tolerance"] = self.tolerance * tolerance_scale
        config_dict["solver"]["convergence"] = suite_criterion

        cfg = pyamgx.Config().create_from_dict(config_dict)
        solver = pyamgx.Solver().create(self.resources, cfg)
        solver.setup(self.A_gpu)

        # solve system on GPU

        self._log.debug("BEGIN solve")

        solver.solve(self.b_gpu, self.x_gpu)

        self._log.debug("END solve")

        # download values from GPU to CPU
        self.x_gpu.download(x)

        if solver.iterations_number == -1:
            residual = None
        else:
            residual = solver.get_residual()

        self._setConvergence(suite="pyamgx",
                             code=solver.status,
                             iterations=solver.iterations_number,
                             tolerance_scale=tolerance_scale,
                             residual=residual)

        self.convergence.warn()

        solver.destroy()
        cfg.destroy()

        return x

    @property
    def _Lxb(self):
        """Matrix, solution vector, and right-hand side vector

        Returns
        -------
        L : ~fipy.matrices.scipyMatrix._ScipyMeshMatrix
            Sparse matrix object
        x : ndarray
            Solution variable
        b : ndarray
            Right-hand side vector
        """
        return (self.matrix, self.var.ravel(), numerix.array(self.RHSvector))

    def _solve(self):
         if self.var.mesh.communicator.Nproc > 1:
             raise Exception("pyamgx solvers cannot be used with multiple processors")

         self.var[:] = numerix.reshape(self._solve_(self.matrix, self.var.ravel(), numerix.array(self.RHSvector)), self.var.shape)
