from __future__ import unicode_literals
import numpy
from scipy.sparse import csr_matrix

import pyamgx

from fipy.solvers.solver import Solver
from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.tools import numerix
from fipy.tools.timer import Timer

__all__ = ["PyAMGXSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PyAMGXSolver(Solver):

    def __init__(self, config_dict, tolerance=1e-10, iterations=2000,
                 precon=None, smoother=None, **kwargs):
        """
        Parameters
        ----------
        config_dict : dict
            AMGX configuration options
        tolerance : float
            Required error tolerance.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pyamgx.preconditioners.preconditioners.Preconditioner, optional
        smoother : ~fipy.solvers.pyamgx.smoothers.smoothers.Smoother, optional
        **kwargs
            Other AMGX solver options
        """
        # update solver config:
        config_dict["solver"]["tolerance"] = tolerance
        config_dict["solver"]["max_iters"] = iterations
        if precon:
            config_dict["solver"]["preconditioner"] = precon
        if smoother:
            config_dict["solver"]["smoother"] = smoother
        config_dict["solver"].update(kwargs)

        # create AMGX objects:
        self.cfg = pyamgx.Config().create_from_dict(config_dict)
        self.resources = pyamgx.Resources().create_simple(self.cfg)
        self.x_gpu = pyamgx.Vector().create(self.resources)
        self.b_gpu = pyamgx.Vector().create(self.resources)
        self.A_gpu = pyamgx.Matrix().create(self.resources)
        self.solver = pyamgx.Solver().create(self.resources, self.cfg)

        super(PyAMGXSolver, self).__init__(tolerance=tolerance, iterations=iterations)

    def __exit__(self, *args):
        # destroy AMGX objects:
        self.A_gpu.destroy()
        self.b_gpu.destroy()
        self.x_gpu.destroy()
        self.solver.destroy()
        self.resources.destroy()
        self.cfg.destroy()

    @property
    def _matrixClass(self):
        return _ScipyMeshMatrix

    def _storeMatrix(self, var, matrix, RHSvector):
        self.var = var
        self.matrix = matrix
        self.RHSvector = RHSvector
        self.A_gpu.upload_CSR(self.matrix.matrix)
        self.solver.setup(self.A_gpu)

    def _solve_(self, L, x, b):
        # transfer data from CPU to GPU
        self._log.debug("BEGIN cpu2gpu")

        with Timer() as t:
            self.x_gpu.upload(x)
            self.b_gpu.upload(b)

        self._log.debug("END cpu2gpu - {elapsed} ns".format(elapsed=t.elapsed))

        # solve system on GPU
        self._log.debug("BEGIN solve")

        with Timer() as t:
            self.solver.solve(self.b_gpu, self.x_gpu)

        self._log.debug("END solve - {} ns".format(t.elapsed))

        # download values from GPU to CPU
        self._log.debug("BEGIN gpu2cpu")

        with Timer() as t:
            self.x_gpu.download(x)

        self._log.debug("END gpu2cpu - {} ns".format(t.elapsed))

        return x

    def _solve(self):
         if self.var.mesh.communicator.Nproc > 1:
             raise Exception("SciPy solvers cannot be used with multiple processors")

         self.var[:] = numerix.reshape(self._solve_(self.matrix, self.var.ravel(), numerix.array(self.RHSvector)), self.var.shape)
