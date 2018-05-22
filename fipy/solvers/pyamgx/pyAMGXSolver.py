#!/usr/bin/env python
import numpy
from scipy.sparse import csr_matrix

import pyamgx

from fipy.solvers.solver import Solver
from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.tools import numerix

__all__ = ["PyAMGXSolver"]

class PyAMGXSolver(Solver):

    def __init__(self, config_dict, tolerance=1e-10, iterations=2000,
                 preconditioner=None, smoother=None, **kwargs):
        """
        :Parameters:
          - `config_dict`: Dictionary specifying AMGX configuration options.
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `preconditioner`: Preconditioner to use.
          - `smoother`: Smoother to use.
          - `kwargs` - Keyword arguments specifying other solver options (see AMGX reference guide).
        """
        # update solver config:
        config_dict["solver"]["tolerance"] = tolerance
        config_dict["solver"]["max_iters"] = iterations
        if preconditioner:
            config_dict["solver"]["preconditioner"] = preconditioner
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
        self.x_gpu.upload(x)
        self.b_gpu.upload(b)

        # solve system on GPU
        self.solver.solve(self.b_gpu, self.x_gpu)

        # download values from GPU to CPU
        self.x_gpu.download(x)
        return x

    def _solve(self):
         if self.var.mesh.communicator.Nproc > 1:
             raise Exception("SciPy solvers cannot be used with multiple processors")

         self.var[:] = numerix.reshape(self._solve_(self.matrix, self.var.ravel(), numerix.array(self.RHSvector)), self.var.shape)
