
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
