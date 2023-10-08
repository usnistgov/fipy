from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools.numerix import random
from fipy.variables.noiseVariable import NoiseVariable

__all__ = ["UniformNoiseVariable"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class UniformNoiseVariable(NoiseVariable):
    r"""
    Represents a uniform distribution of random numbers.

    We generate noise on a uniform Cartesian mesh

    >>> from fipy.meshes import Grid2D
    >>> noise = UniformNoiseVariable(mesh=Grid2D(nx=100, ny=100))

    and histogram the noise

    >>> from fipy.variables.histogramVariable import HistogramVariable
    >>> histogram = HistogramVariable(distribution=noise, dx=0.01, nx=120, offset=-.1)

    >>> if __name__ == '__main__':
    ...     from fipy import Viewer
    ...     viewer = Viewer(vars=noise,
    ...                     datamin=0, datamax=1)
    ...     histoplot = Viewer(vars=histogram)

    >>> from builtins import range
    >>> for i in range(10):
    ...     noise.scramble()
    ...     if __name__ == '__main__':
    ...         viewer.plot()
    ...         histoplot.plot()

    .. image:: /figures/fipy/variables/uniform.*
       :scale: 25
       :align: center
       :alt: random values with a uniform distribution

    .. image:: /figures/fipy/variables/uni-histogram.*
       :scale: 25
       :align: center
       :alt: histogram of random values with a uniform distribution
    """
    def __init__(self, mesh, name = '', minimum = 0., maximum = 1., hasOld = 0):
        """
        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The mesh on which to define the noise.
        minimum : float
            The minimum (not-inclusive) value of the distribution.
        maximum : float
            The maximum (not-inclusive) value of the distribution.
        """
        self.minimum = minimum
        self.maximum = maximum
        NoiseVariable.__init__(self, mesh = mesh, name = name, hasOld = hasOld)

    def random(self):
        return random.uniform(self.minimum, self.maximum,
                              size=[self.mesh.globalNumberOfCells])

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
