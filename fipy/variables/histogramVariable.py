## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "histogramVariable.py"
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
