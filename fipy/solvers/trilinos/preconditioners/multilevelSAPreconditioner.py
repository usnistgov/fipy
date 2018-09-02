#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "multilevelSAPreconditioner.py"
 #                                    created: 06/25/07
 #                                last update: 11/8/11 {4:09:31 PM}
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Maxsim Gibiansky <maxsim.gibiansky@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 #  Description:
 #
 #  History
 #
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2007-06-25 MLG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelSAPreconditioner"]

class MultilevelSAPreconditioner(Preconditioner):
    """
    Multilevel preconditioner for Trilinos solvers suitable classical
    smoothed aggregation for symmetric positive definite or nearly
    symmetric positive definite systems.
    """

    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList({"output": 0,
                                    "max levels" : 10,
                                    "prec type" : "MGV",
                                    "increasing or decreasing" : "increasing",
                                    "aggregation: type" : "Uncoupled-MIS",
                                    "aggregation: damping factor" : 4. / 3.,
##                                    "energy minimization: enable" : False,
##                                    "smoother: type" : "Aztec",
##                                    "smoother: type" : "symmetric Gauss-Seidel",
##                                    "eigen-analysis: type" : "power-method",
                                    "eigen-analysis: type" : "cg",
                                    "eigen-analysis: iterations" : 10,
                                    "smoother: sweeps" : 2,
                                    "smoother: damping factor" : 1.0,
                                    "smoother: pre or post" : 'both',
                                    "smoother: type" : "symmetric Gauss-Seidel",
                                    "coarse: type" : 'Amesos-KLU',
                                    "coarse: max size" : 128
                                    })

        self.Prec.ComputePreconditioner()

        solver.SetPrecOperator(self.Prec)
