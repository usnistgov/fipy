from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.viewers.matplotlibViewer.abstractMatplotlibViewer import AbstractMatplotlibViewer

class AbstractMatplotlib2DViewer(AbstractMatplotlibViewer):
    """Base class for plotting 2D :class:`~fipy.variables.meshVariable.MeshVariable` objects with Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    def figaspect(self, figaspect, colorbar):
        if figaspect == 'auto':
            figaspect = self.vars[0].mesh.aspect2D
            # We can't make the colorbar, yet (as we don't even have a figure)
            # so hardcode these values, based on
            # https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html#matplotlib.figure.Figure.colorbar
            #   fraction	0.15; fraction of original axes to use for colorbar
            #   pad	        0.05 if vertical, 0.15 if horizontal; fraction of original axes between colorbar and new image axes
            if colorbar == "vertical":
                figaspect = figaspect * (1. - (0.15 + 0.05))
            elif colorbar == "horizontal":
                figaspect = figaspect / (1. - (0.15 + 0.15))
        return figaspect

    def _plot(self):
        zmin, zmax = self._autoscale(vars=self.vars,
                                     datamin=self._getLimit(('datamin', 'zmin')),
                                     datamax=self._getLimit(('datamax', 'zmax')))

        self._norm.vmin = zmin
        self._norm.vmax = zmax

        self._mappable.set_norm(self._norm)
