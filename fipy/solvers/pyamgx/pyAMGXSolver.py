#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix

import pyamgx

from fipy.solvers.solver import Solver
from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.tools import numerix

__all__ = ["PyAMGXSolver"]

class PyAMGXSolver(Solver):
    """
    The `PyAMGXSolver` is a general interface to pyamgx.
    """
    def __init__(self, config_dict, *args, **kwargs):
        """
        Parameters
        ----------
        config_dict : dict
            Dictionary specifying AMGX configuration options
        """
        self.config_dict = config_dict
        self.cfg = pyamgx.Config().create_from_dict(self.config_dict)
        self.resources = pyamgx.Resources().create_simple(self.cfg)
        self.x_gpu = pyamgx.Vector().create(self.resources)
        self.b_gpu = pyamgx.Vector().create(self.resources)
        self.A_gpu = pyamgx.Matrix().create(self.resources)
        self.solver = pyamgx.Solver().create(self.resources, self.cfg)

    def __exit__(self, exc_type, exc_value, traceback):
        self.solver.destroy()
        self.A_gpu.destroy()
        self.b_gpu.destroy()
        self.x_gpu.destroy()
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

    def _solve(self):
        self._solve_(self.matrix, self.var.ravel(), numerix.array(self.RHSvector))
            
    def _canSolveAsymmetric(self):
        return False
