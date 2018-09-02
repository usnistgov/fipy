#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "commWrapper.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

from fipy.tools import numerix

__all__ = ["CommWrapper", "ParallelCommWrapper"]

class CommWrapper(object):
    """MPI Communicator wrapper

    Encapsulates capabilities needed for Epetra. Some capabilities are not parallel.

    """

    def __init__(self, Epetra=None):
        self.epetra_comm = Epetra.PyComm()

    def __repr__(self):
        return "%s()" % self.__class__.__name__

    @property
    def procID(self):
        return self.epetra_comm.MyPID()

    @property
    def Nproc(self):
        return self.epetra_comm.NumProc()

    def Barrier(self):
        self.epetra_comm.Barrier()

    def all(self, a, axis=None):
        return a.all(axis=axis)

    def any(self, a, axis=None):
        return a.any(axis=axis)

    def allclose(self, a, b, rtol=1.e-5, atol=1.e-8):
        return numerix.allclose(a, b, rtol=rtol, atol=atol)

    def allequal(self, a, b):
        return numerix.allequal(a, b)

    def bcast(self, obj, root=0):
        return obj

    def allgather(self, obj):
        return obj

    def sum(self, a, axis=None):
        summed = numerix.array(a).sum(axis=axis)
        shape = summed.shape
        if shape == ():
            summed = summed.reshape((1,))
        parallelSummed = self.epetra_comm.SumAll(summed)
        if shape == ():
            parallelSummed = parallelSummed.reshape(())
        return parallelSummed

    def __getstate__(self):
        return {'dummy': 0}

    def __setstate__(self, dict):
        from PyTrilinos import Epetra
        self.__init__(Epetra=Epetra)

    def Norm2(self, vec):
        return vec.Norm2()

    def MaxAll(self, vec):
        return self.epetra_comm.MaxAll(numerix.array(vec))

    def MinAll(self, vec):
        return self.epetra_comm.MinAll(numerix.array(vec))

class ParallelCommWrapper(CommWrapper):
    """MPI Communicator wrapper for parallel processes"""
    pass
