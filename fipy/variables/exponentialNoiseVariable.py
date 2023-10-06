from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

from fipy.tools.numerix import random
from fipy.variables.noiseVariable import NoiseVariable

__all__ = ["ExponentialNoiseVariable"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ExponentialNoiseVariable(NoiseVariable):
    r"""
    Represents an exponential distribution of random numbers with the probability
    distribution

    .. math::

       \mu^{-1} e^{-\frac{x}{\mu}}

    with a mean parameter :math:`\mu`.

    Seed the random module for the sake of deterministic test results.

    >>> from fipy import numerix
    >>> numerix.random.seed(1)

    We generate noise on a uniform Cartesian mesh

    >>> from fipy.variables.variable import Variable
    >>> mean = Variable()
    >>> from fipy.meshes import Grid2D
    >>> noise = ExponentialNoiseVariable(mesh = Grid2D(nx = 100, ny = 100), mean = mean)

    We histogram the root-volume-weighted noise distribution

    >>> from fipy.variables.histogramVariable import HistogramVariable
    >>> histogram = HistogramVariable(distribution = noise, dx = 0.1, nx = 100)

    and compare to a Gaussian distribution

    >>> from fipy.variables.cellVariable import CellVariable
    >>> expdist = CellVariable(mesh = histogram.mesh)
    >>> x = histogram.mesh.cellCenters[0]

    >>> if __name__ == '__main__':
    ...     from fipy import Viewer
    ...     viewer = Viewer(vars=noise, datamin=0, datamax=5)
    ...     histoplot = Viewer(vars=(histogram, expdist),
    ...                        datamin=0, datamax=1.5)

    >>> from fipy.tools.numerix import arange, exp

    >>> for mu in arange(0.5, 3, 0.5):
    ...     mean.value = (mu)
    ...     expdist.value = ((1/mean)*exp(-x/mean))
    ...     if __name__ == '__main__':
    ...         import sys
    ...         print("mean: %g" % mean, file=sys.stderr)
    ...         viewer.plot()
    ...         histoplot.plot()

    >>> print(abs(noise.faceGrad.divergence.cellVolumeAverage) < 5e-15)
    1

    .. image:: /figures/fipy/variables/exp.*
      :scale: 25
      :align: center
      :alt: random values with an exponential distribution

    .. image:: /figures/fipy/variables/exp-histogram.*
      :scale: 25
      :align: center
      :alt: histogram of random values with an exponential distribution

    """
    def __init__(self, mesh, mean=0.0, name = '', hasOld = 0):
        r"""
        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The mesh on which to define the noise.
        mean : float
            The mean of the distribution :math:`\mu`.
        """
        NoiseVariable.__init__(self, mesh = mesh, name = name, hasOld = hasOld)
        self.mean = self._requires(mean)

    def random(self):
        return random.exponential(scale = self.mean,
                                  size = [self.mesh.globalNumberOfCells])

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


