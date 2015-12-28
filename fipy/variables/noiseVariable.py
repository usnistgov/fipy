## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "noiseVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  noiseVariable.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 #
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution of
 #  this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable

__all__ = ["NoiseVariable"]

class NoiseVariable(CellVariable):
    r"""
    .. attention:: This class is abstract. Always create one of its subclasses.

    A generic base class for sources of noise distributed over the cells of a mesh.

    In the event that the noise should be conserved, use::

        <Specific>NoiseVariable(...).faceGrad.divergence

    The `seed()` and `get_seed()` functions of the
    `fipy.tools.numerix.random` module can be set and query the random
    number generated used by all `NoiseVariable` objects.
    """
    def __init__(self, mesh, name = '', hasOld = 0):
        if self.__class__ is NoiseVariable:
            raise NotImplementedError, "can't instantiate abstract base class"

        CellVariable.__init__(self, mesh = mesh, name = name, hasOld = hasOld)
        self.scramble()

    def copy(self):
        """
        Copy the value of the `NoiseVariable` to a static `CellVariable`.
        """
        return CellVariable(mesh = self.mesh,
                            name = self.name + "_old",
                            value = self.value,
                            hasOld = 0)

    def scramble(self):
        """
        Generate a new random distribution.
        """
        self._markStale()

    def random(self):
        pass

    def parallelRandom(self):

        if self.mesh.communicator.procID == 0:
            return self.random()
        else:
            return None

    def _calcValue(self):
        from fipy.tools import parallelComm

        rnd = self.parallelRandom()

        if parallelComm.Nproc > 1:
            rnd = parallelComm.bcast(rnd, root=0)

            return rnd[self.mesh._globalOverlappingCellIDs]
        else:
            return rnd
