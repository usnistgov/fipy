__docformat__ = 'restructuredtext'

__all__ = ["AbstractMatplotlibViewer"]

from fipy.viewers.viewer import AbstractViewer

def _isnotebook():
    """return True if running in a jupyter notebook
     https://stackoverflow.com/a/39662359/2019542
    """
    try:
        import IPython

        shell = IPython.get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except (ImportError, NameError):
        return False      # Probably standard Python interpreter

def _isretina():
    """return True if jupyter notebook configuration
    InlineBackend.figure_formats contains 'retina' or
    InlineBackend.figure_format is set to 'retina'
    """
    isretina = False

    try:
        import IPython

        cfg = IPython.get_ipython().config
        if "InlineBackend" in cfg:
            inbk = cfg["InlineBackend"]
            if "figure_formats" in inbk:
                if "retina" in inbk["figure_formats"]:
                    isretina = True
            if not isretina and ("figure_format" in inbk):
                isretina = inbk["figure_format"] == "retina"
    except ImportError:
        pass
    return isretina

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

        from matplotlib import pyplot as plt

        plt.ion()

        if axes is None:
            w, h = plt.figaspect(self.figaspect(figaspect))
            fig = plt.figure(figsize=(w, h))
            self.axes = plt.gca()
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
        from matplotlib import pyplot as plt

        fig = self.axes.get_figure()

        plt.ioff()

        self._plot()

        plt.draw()

        try:
            fig.canvas.flush_events()
        except NotImplementedError:
            pass

        plt.ion()

        if _isnotebook():
            # plots don't animate in the notebook unless we
            # explicitly clear_output and display
            from IPython.display import display, clear_output

            clear_output(wait=True)
            display(self)
        else:
            fig.show()

        if filename is not None:
            fig.savefig(filename)

    def _validFileExtensions(self):
        return ["""
        Matplotlib has no reliable way to determine
        valid file extensions. Either guess, or see
        <http://matplotlib.sourceforge.net/faq/installing_faq.html#backends>
        and then guess. Yes, this is lame.
        """]

#         filetypes = plt.figure(self.id).canvas.filetypes
#         return [".%s" % key for key in filetypes.keys()]

    def _repr_png_(self):
        """Render as a PNG for IPython notebook, per display_protocol.ipynb

        Invoke with `display(myViewer)`
        """
        from IPython.core.pylabtools import print_figure, retina_figure

        fig = self.axes.get_figure()
        if _isretina():
            return retina_figure(fig)
        else:
            return print_figure(fig, "png")

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
