#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "betaNoiseVariable.py"
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

__all__ = ["BetaNoiseVariable"]

class BetaNoiseVariable(NoiseVariable):
    r"""
    Represents a beta distribution of random numbers with the probability
    distribution

    .. math::

       x^{\alpha - 1}\frac{\beta^\alpha e^{-\beta x}}{\Gamma(\alpha)}

    with a shape parameter :math:`\alpha`, a rate parameter :math:`\beta`, and
    :math:`\Gamma(z) = \int_0^\infty t^{z - 1}e^{-t}\,dt`.

    Seed the random module for the sake of deterministic test results.

    >>> from fipy import numerix
    >>> numerix.random.seed(1)

    We generate noise on a uniform Cartesian mesh

    >>> from fipy.variables.variable import Variable
    >>> alpha = Variable()
    >>> beta = Variable()
    >>> from fipy.meshes import Grid2D
    >>> noise = BetaNoiseVariable(mesh = Grid2D(nx = 100, ny = 100), alpha = alpha, beta = beta)

    We histogram the root-volume-weighted noise distribution

    >>> from fipy.variables.histogramVariable import HistogramVariable
    >>> histogram = HistogramVariable(distribution = noise, dx = 0.01, nx = 100)

    and compare to a Gaussian distribution

    >>> from fipy.variables.cellVariable import CellVariable
    >>> betadist = CellVariable(mesh = histogram.mesh)
    >>> x = CellVariable(mesh=histogram.mesh, value=histogram.mesh.cellCenters[0])
    >>> from scipy.special import gamma as Gamma # doctest: +SCIPY
    >>> betadist = ((Gamma(alpha + beta) / (Gamma(alpha) * Gamma(beta)))
    ...             * x**(alpha - 1) * (1 - x)**(beta - 1)) # doctest: +SCIPY

    >>> if __name__ == '__main__':
    ...     from fipy import Viewer
    ...     viewer = Viewer(vars=noise, datamin=0, datamax=1)
    ...     histoplot = Viewer(vars=(histogram, betadist),
    ...                        datamin=0, datamax=1.5)

    >>> from fipy.tools.numerix import arange

    >>> for a in arange(0.5,5,0.5):
    ...     alpha.value = a
    ...     for b in arange(0.5,5,0.5):
    ...         beta.value = b
    ...         if __name__ == '__main__':
    ...             import sys
    ...             print >>sys.stderr, "alpha: %g, beta: %g" % (alpha, beta)
    ...             viewer.plot()
    ...             histoplot.plot()

    >>> print abs(noise.faceGrad.divergence.cellVolumeAverage) < 5e-15
    1

    .. image:: fipy/variables/beta.*
      :scale: 25
      :align: center
      :alt: random values with a beta distribution

    .. image:: fipy/variables/beta-histogram.*
      :scale: 25
      :align: center
      :alt: histogram of random values with a beta distribution

    """
    def __init__(self, mesh, alpha, beta, name = '', hasOld = 0):
        r"""
        :Parameters:
            - `mesh`: The mesh on which to define the noise.
            - `alpha`: The parameter :math:`\alpha`.
            - `beta`: The parameter :math:`\beta`.

        """
        NoiseVariable.__init__(self, mesh = mesh, name = name, hasOld = hasOld)
        self.alpha = self._requires(alpha)
        self.beta = self._requires(beta)

    def random(self):
        return random.beta(a = self.alpha, b = self.beta,
                           size = [self.mesh.globalNumberOfCells])

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
