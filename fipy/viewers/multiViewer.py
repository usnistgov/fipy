## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "multiViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
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

from fipy.viewers.viewer import AbstractViewer

__all__ = ["MultiViewer"]

class MultiViewer(AbstractViewer):
    """
    Treat a collection of different viewers (such for different 2D plots
    or 1D plots with different axes) as a single viewer that will `plot()`
    all subviewers simultaneously.
    """
    def __init__(self, viewers):
        """
        :Parameters:
          viewers : list
            the viewers to bind together
        """
        if type(viewers) not in [type([]), type(())]:
            viewers = [viewers]
        self.viewers = viewers

    def setLimits(self, limits={}, **kwlimits):
        kwlimits.update(limits)
        for viewer in self.viewers:
            viewer.setLimits(**kwlimits)

    def plot(self):
        for viewer in self.viewers:
            viewer.plot()
