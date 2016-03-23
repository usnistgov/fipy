#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "matplotlibViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = ["AbstractMatplotlibViewer"]

from fipy.viewers.viewer import AbstractViewer

class AbstractMatplotlibViewer(AbstractViewer):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    The `AbstractMatplotlibViewer` is the base class for the viewers that use the
    Matplotlib_ python plotting package.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    def __init__(self, vars, title=None, figaspect=1.0, cmap=None, colorbar=None, axes=None, log=False, **kwlimits):
        """
        Create a `AbstractMatplotlibViewer`.

        :Parameters:
          vars
            a `CellVariable` or tuple of `CellVariable` objects to plot
          title
            displayed at the top of the `Viewer` window
          figaspect
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is an array, figaspect will determine the width and
            height for a figure that would fit array preserving aspect ratio.
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.
          cmap
            the colormap. Defaults to `matplotlib.cm.jet`
          colorbar
            plot a colorbar in specified orientation if not `None`
          axes
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
          log
            whether to logarithmically scale the data
        """
        if self.__class__ is AbstractMatplotlibViewer:
            raise NotImplementedError("can't instantiate abstract base class")

        AbstractViewer.__init__(self, vars=vars, title=title, **kwlimits)

        import pylab

        pylab.ion()

        if axes is None:
            w, h = pylab.figaspect(self.figaspect(figaspect))
            fig = pylab.figure(figsize=(w, h))
            self.axes = pylab.gca()
        else:
            self.axes = axes
            fig = axes.get_figure()

        self.id = fig.number

        self.axes.set_title(self.title)

        import matplotlib
        # Set the colormap and norm to correspond to the data for which
        # the colorbar will be used.
        if cmap is None:
            self.cmap = matplotlib.cm.jet
        else:
            self.cmap = cmap

        if colorbar:
            self.colorbar = _ColorBar(viewer=self)
        else:
            self.colorbar = None

        self.norm = None
        self.log = log

        try:
            # Plotting needs to work differently for inline
            # integration in the IPython notebook.
            # (test is from http://stackoverflow.com/a/15346737/2019542)
            backend = pylab.get_backend()
            self.IPYinline = __IPYTHON__ and ("inline" in backend)
        except NameError:
            self.IPYinline = False

    def figaspect(self, figaspect):
        return figaspect

    def log():
        doc = "logarithmic data scaling"

        def fget(self):
            from matplotlib import colors
            return isinstance(self.norm, colors.LogNorm)

        def fset(self, value):
            from matplotlib import colors
            if value:
                self.norm = colors.LogNorm()
            else:
                self.norm = colors.Normalize()

            if self.colorbar is not None:
                self.colorbar.set_norm(self.norm)

        return locals()

    log = property(**log())

    def plot(self, filename = None):
        import pylab

        fig = pylab.figure(self.id)

        if self.IPYinline:
            from IPython.display import clear_output, display_png

            clear_output(wait=True)
            display_png(self)
        else:
            pylab.ioff()

            self._plot()

            pylab.draw()

            try:
                fig.canvas.flush_events()
            except NotImplementedError:
                pass

            pylab.ion()

            pylab.show(block=False)

        if filename is not None:
            pylab.savefig(filename)

    def _validFileExtensions(self):
        import pylab
        return ["""
        Matplotlib has no reliable way to determine
        valid file extensions. Either guess, or see
        <http://matplotlib.sourceforge.net/faq/installing_faq.html#backends>
        and then guess. Yes, this is lame.
        """]

#         filetypes = pylab.figure(self.id).canvas.filetypes
#         return [".%s" % key for key in filetypes.keys()]

    def _repr_png_(self):
        """Render as a PNG for IPython notebook, per display_protocol.ipynb

        Invoke with `display(myViewer)`
        """
        from IPython.core.pylabtools import print_figure

        self._plot()
        return print_figure(fig=self.axes.get_figure(), fmt="png")

class _ColorBar(object):
    def __init__(self, viewer, vmin=-1, vmax=1, orientation="vertical"):
        self.viewer = viewer

        import matplotlib
        cbax, kw = matplotlib.colorbar.make_axes(viewer.axes, orientation=orientation)

        # ColorbarBase derives from ScalarMappable and puts a colorbar
        # in a specified axes, so it has everything needed for a
        # standalone colorbar.  There are many more kwargs, but the
        # following gives a basic continuous colorbar with ticks
        # and labels.
        import matplotlib.colors as colors
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        self._cb = matplotlib.colorbar.ColorbarBase(cbax, norm=norm, cmap=viewer.cmap,
                                                    orientation=orientation)
        self._cb.set_label(viewer.vars[0].name)

        self.formatter = None

    def get_norm(self):
        return self._cb.get_norm()

    def set_norm(self, value):
        self._cb.set_norm(value)
        if self.formatter is None:
            from matplotlib import colors, ticker
            if isinstance(value, colors.LogNorm):
                self._cb.formatter = ticker.LogFormatterMathtext()
            else:
                self._cb.formatter = ticker.ScalarFormatter()

    norm = property(fget=get_norm, fset=set_norm, doc="data normalization")

    def plot(self): #, vmin, vmax):
        self._cb.set_norm(self.viewer.norm)
        self._cb.cmap = self.viewer.cmap
        self._cb.draw_all()

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
