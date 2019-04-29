from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable

__all__ = ["NoiseVariable"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

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
            raise NotImplementedError("can't instantiate abstract base class")

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
