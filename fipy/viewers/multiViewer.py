from __future__ import unicode_literals
from fipy.viewers.viewer import AbstractViewer

__all__ = ["MultiViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MultiViewer(AbstractViewer):
    """
    Treat a collection of different viewers (such for different 2D plots
    or 1D plots with different axes) as a single viewer that will `plot()`
    all subviewers simultaneously.
    """
    def __init__(self, viewers):
        """
        Parameters
        ----------
        viewers : :obj:`list` of ~fipy.viewers.viewer.Viewer
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
