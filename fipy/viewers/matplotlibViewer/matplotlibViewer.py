from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = ["AbstractMatplotlibViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

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
    """return `True` if jupyter notebook configuration
    `InlineBackend.figure_formats` contains `retina` or
    `InlineBackend.figure_format` is set to `retina`
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

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or list
            the `Variable` objects to display.
        title : str, optional
            displayed at the top of the `Viewer` window
        figaspect : float, optional
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is `auto`, the aspect ratio will be determined from
            the Variable's mesh.
        float xmin, xmax, ymin, ymax, datamin, datamax : float, optional
            displayed range of data. A 1D `Viewer` will only use *xmin* and
            *xmax*, a 2D viewer will also use *ymin* and *ymax*. All
            viewers will use *datamin* and *datamax*. Any limit set to a
            (default) value of `None` will autoscale.
        cmap : ~matplotlib.colors.Colormap, optional
            the :class:`~matplotlib.colors.Colormap`.
            Defaults to `matplotlib.cm.jet`
        colorbar : bool, optional
            plot a color bar in specified orientation if not `None`
        axes : ~matplotlib.axes.Axes, optional
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        log : bool, optional
          whether to logarithmically scale the data
        """
        if self.__class__ is AbstractMatplotlibViewer:
            raise NotImplementedError("can't instantiate abstract base class")

        AbstractViewer.__init__(self, vars=vars, title=title, **kwlimits)

        from matplotlib import pyplot as plt

        plt.ion()

        if axes is None:
            w, h = plt.figaspect(self.figaspect(figaspect, colorbar))
            self.fig = plt.figure(figsize=(w, h))
            self.axes = plt.gca()
        else:
            self.axes = axes
            self.fig = axes.get_figure()

        self.id = self.fig.number

        self.axes.set_title(self.title)

        self._mappable = None

        import matplotlib
        # Set the colormap and norm to correspond to the data for which
        # the colorbar will be used.
        if cmap is None:
            self.cmap = matplotlib.cm.jet
        else:
            self.cmap = cmap

        self.norm = None
        self.log = log

        if colorbar:
            self.colorbar = self.fig.colorbar(mappable=self.mappable,
                                              orientation=colorbar,
                                              label=self.vars[0].name)
        else:
            self.colorbar = None

    def figaspect(self, figaspect, colorbar):
        return figaspect

    def _make_mappable(self):
        import matplotlib
        mappable = matplotlib.cm.ScalarMappable(norm=self.norm, cmap=self.cmap)
        # ignored, but needed for matplotlib < 3.0 (Py2k)
        mappable.set_array(self.vars[0].value)

        return mappable

    @property
    def mappable(self):
        if self._mappable is None:
            self._mappable = self._make_mappable()
        return self._mappable

    def log():
        doc = "logarithmic data scaling"

        def fget(self):
            from matplotlib import colors
            return isinstance(self.norm, colors.LogNorm)

        def fset(self, value):
            zmin, zmax = self._autoscale(vars=self.vars,
                                         datamin=self._getLimit(('datamin', 'zmin')),
                                         datamax=self._getLimit(('datamax', 'zmax')))

            from matplotlib import colors
            if value:
                self.norm = colors.LogNorm(vmin=zmin, vmax=zmax)
            else:
                self.norm = colors.Normalize(vmin=zmin, vmax=zmax)

            self.mappable.set_norm(self.norm)

        return locals()

    log = property(**log())

    def plot(self, filename = None):
        from matplotlib import pyplot as plt

        plt.ioff()

        self._plot()

        plt.draw()

        try:
            self.fig.canvas.flush_events()
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
            self.fig.show()

        if filename is not None:
            self.fig.savefig(filename)

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
        """Render as a PNG for IPython notebook, per `display_protocol.ipynb`

        Invoke with `display(myViewer)`
        """
        from IPython.core.pylabtools import print_figure, retina_figure

        if _isretina():
            return retina_figure(self.fig)
        else:
            return print_figure(self.fig, "png")

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
