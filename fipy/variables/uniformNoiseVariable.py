## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "uniformNoiseVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed by employees of the National Institute
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # works of NIST employees are not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
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

from fipy.tools.numerix import random
from fipy.variables.noiseVariable import NoiseVariable

__all__ = ["UniformNoiseVariable"]

class UniformNoiseVariable(NoiseVariable):
    r"""
    Represents a uniform distribution of random numbers.

    We generate noise on a uniform cartesian mesh

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

    >>> for i in range(10):
    ...     noise.scramble()
    ...     if __name__ == '__main__':
    ...         viewer.plot()
    ...         histoplot.plot()

    .. image:: fipy/variables/uniform.*
       :scale: 25
       :align: center
       :alt: random values with a uniform distribution

    .. image:: fipy/variables/uni-histogram.*
       :scale: 25
       :align: center
       :alt: histogram of random values with a uniform distribution
    """
    def __init__(self, mesh, name = '', minimum = 0., maximum = 1., hasOld = 0):
        """
        :Parameters:
            - `mesh`: The mesh on which to define the noise.
            - `minimum`: The minimum (not-inclusive) value of the distribution.
            - `maximum`: The maximum (not-inclusive) value of the distribution.
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
