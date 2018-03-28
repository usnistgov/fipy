#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "multilevelDDPreconditioner.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelDDPreconditioner"]

class MultilevelDDPreconditioner(Preconditioner):
    """
    Multilevel preconditioner for Trilinos solvers. A classical smoothed
    aggregation-based 2-level domain decomposition.
    """

    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList({"output": 0,
                                    "max levels" : 2,
                                    "prec type" : "MGV",
                                    "increasing or decreasing" : "increasing",
                                    "aggregation: type" : "METIS",
                                    "aggregation: local aggregates" : 1,
                                    "aggregation: damping factor" : 4. / 3.,
                                    "eigen-analysis: type" : "power-method",
                                    "eigen-analysis: iterations" : 20,
                                    "smoother: sweeps" : 1,
                                    "smoother: pre or post" : 'both',
                                    "smoother: type" : "Aztec",
                                    "smoother: Aztec as solver" : False,
                                    "coarse: type" : 'Amesos-KLU',
                                    "coarse: max size" : 128
                                    })

        self.Prec.ComputePreconditioner()

        solver.SetPrecOperator(self.Prec)
