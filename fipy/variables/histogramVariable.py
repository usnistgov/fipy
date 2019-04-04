__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.meshes import Grid1D
from fipy.tools import numerix

__all__ = ["HistogramVariable"]

class HistogramVariable(CellVariable):
    def __init__(self, distribution, dx = 1., nx = None, offset = 0.):
        r"""
        Produces a histogram of the values of the supplied distribution.

        :Parameters:

            - `distribution`: The collection of values to sample.
            - `dx`: the bin size
            - `nx`: the number of bins
            - `offset`: the position of the first bin
        """
        CellVariable.__init__(self, mesh = Grid1D(dx = dx, nx = nx) + (offset,))
        self.distribution = self._requires(distribution)

    def _calcValue(self):
        l = len(self.distribution)
        bins = self.mesh.cellCenters[0]
        n = numerix.searchsorted(numerix.sort(self.distribution), bins)
        n = numerix.concatenate([n, [l]])
        dx = bins[1:] - bins[:-1]
        return (n[1:] - n[:-1]) / numerix.concatenate([dx, [dx[-1]]]) / float(l)
