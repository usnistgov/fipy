from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools.numerix import random, sqrt
from fipy.variables.noiseVariable import NoiseVariable

__all__ = ["GaussianNoiseVariable"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class GaussianNoiseVariable(NoiseVariable):
    r"""

    Represents a normal (Gaussian) distribution of random numbers with
    mean :math:`\mu` and variance
    :math:`\langle \eta(\vec{r}) \eta(\vec{r}\,') \rangle = \sigma^2`,
    which has the probability distribution

    .. math::

       \frac{1}{\sigma\sqrt{2\pi}} \exp -\frac{(x - \mu)^2}{2\sigma^2}

    For example, the variance of thermal noise that is uncorrelated in space and
    time is often expressed as

    .. math::

       \left\langle
           \eta(\vec{r}, t) \eta(\vec{r}\,', t')
       \right\rangle =
       M k_B T \delta(\vec{r} - \vec{r}\,')\delta(t - t')

    which can be obtained with::

        sigmaSqrd = Mobility * kBoltzmann * Temperature / (mesh.cellVolumes * timeStep)
        GaussianNoiseVariable(mesh = mesh, variance = sigmaSqrd)

    .. note::

       If the time step will change as the simulation progresses, either through
       use of an adaptive iterator or by making manual changes at different
       stages, remember to declare `timeStep` as a `Variable` and to change its
       value with its `setValue()` method.

    >>> import sys
    >>> from fipy.tools.numerix import *

    >>> mean = 0.
    >>> variance = 4.

    Seed the random module for the sake of deterministic test results.

    >>> from fipy import numerix
    >>> numerix.random.seed(3)

    We generate noise on a non-uniform Cartesian mesh with cell dimensions of
    :math:`x^2` and :math:`y^3`.

    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = arange(0.1, 5., 0.1)**2, dy = arange(0.1, 3., 0.1)**3)
    >>> from fipy.variables.cellVariable import CellVariable
    >>> volumes = CellVariable(mesh=mesh, value=mesh.cellVolumes)
    >>> noise = GaussianNoiseVariable(mesh = mesh, mean = mean,
    ...                               variance = variance / volumes)

    We histogram the root-volume-weighted noise distribution

    >>> from fipy.variables.histogramVariable import HistogramVariable
    >>> histogram = HistogramVariable(distribution = noise * sqrt(volumes),
    ...                               dx = 0.1, nx = 600, offset = -30)

    and compare to a Gaussian distribution

    >>> gauss = CellVariable(mesh = histogram.mesh)
    >>> x = histogram.mesh.cellCenters[0]
    >>> gauss.value = ((1/(sqrt(variance * 2 * pi))) * exp(-(x - mean)**2 / (2 * variance)))

    >>> if __name__ == '__main__':
    ...     from fipy.viewers import Viewer
    ...     viewer = Viewer(vars=noise,
    ...                     datamin=-5, datamax=5)
    ...     histoplot = Viewer(vars=(histogram, gauss))

    >>> from builtins import range
    >>> for i in range(10):
    ...     noise.scramble()
    ...     if __name__ == '__main__':
    ...         viewer.plot()
    ...         histoplot.plot()

    >>> print(abs(noise.faceGrad.divergence.cellVolumeAverage) < 5e-15)
    1

    Note that the noise exhibits larger amplitude in the small cells than in the large ones

    .. image:: /figures/fipy/variables/gaussian.*
      :scale: 25
      :align: center
      :alt: random values with a Gaussian distribution

    but that the root-volume-weighted histogram is Gaussian.

    .. image:: /figures/fipy/variables/gauss-histogram.*
      :scale: 25
      :align: center
      :alt: histogram of random values with a Gaussian distribution

    """
    def __init__(self, mesh, name = '', mean = 0., variance = 1., hasOld = 0):
        """
        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The mesh on which to define the noise.
        mean : float
            The mean of the noise distribution, :math:`\mu`.
        variance : float
            The variance of the noise distribution, :math:`\sigma^2`.
        """
        self.mean = mean
        self.variance = variance
        NoiseVariable.__init__(self, mesh = mesh, name = name, hasOld = hasOld)

    def parallelRandom(self):

        if hasattr(self.variance, 'globalValue'):
            variance = self.variance.globalValue
        else:
            variance = self.variance

        if self.mesh.communicator.procID == 0:
            return random.normal(self.mean, sqrt(variance),
                                 size = [self.mesh.globalNumberOfCells])
        else:
            return None

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
