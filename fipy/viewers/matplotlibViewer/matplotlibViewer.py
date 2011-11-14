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

__all__ = ["MatplotlibViewer"]

from fipy.viewers.viewer import _Viewer

def MatplotlibViewer(vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, **kwlimits):
    """Generic function for creating a `MatplotlibViewer`. 
    
    The `MatplotlibViewer` factory will search the module tree and return an
    instance of the first `MatplotlibViewer` it finds of the correct dimension
    and rank.
    
    :Parameters:
      vars
        a `CellVariable` or tuple of `CellVariable` objects to plot
      title
        displayed at the top of the `Viewer` window
      limits : dict
        a (deprecated) alternative to limit keyword arguments
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
        
    It is possible to view different `Variable`s against different Matplotlib_ `Axes
    
    >>> from matplotlib import pylab
    >>> from fipy import *

    >>> fig = pylab.figure()

    >>> ax1 = pylab.subplot((221))
    >>> ax2 = pylab.subplot((223))
    >>> ax3 = pylab.subplot((224))

    >>> k = Variable(name="k", value=0.)

    >>> mesh1 = Grid1D(nx=100)
    >>> x, = mesh1.cellCenters
    >>> xVar = CellVariable(mesh=mesh1, name="x", value=x)
    >>> viewer1 = MatplotlibViewer(vars=(numerix.sin(0.1 * k * xVar), numerix.cos(0.1 * k * xVar / numerix.pi)), 
    ...                            limits={'xmin': 10, 'xmax': 90}, 
    ...                            datamin=-0.9, datamax=2.0,
    ...                            title="Grid1D test",
    ...                            axes=ax1,
    ...                            legend=None)
                                
    >>> mesh2 = Grid2D(nx=50, ny=100, dx=0.1, dy=0.01)
    >>> x, y = mesh2.cellCenters
    >>> xyVar = CellVariable(mesh=mesh2, name="x y", value=x * y)
    >>> viewer2 = MatplotlibViewer(vars=numerix.sin(k * xyVar), 
    ...                            limits={'ymin': 0.1, 'ymax': 0.9}, 
    ...                            datamin=-0.9, datamax=2.0,
    ...                            title="Grid2D test",
    ...                            axes=ax2,
    ...                            colorbar=None)

    >>> mesh3 = (Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
    ...          + (Tri2D(nx=5, ny=5, dx=0.1, dy=0.1) 
    ...             + ((0.5,), (0.2,))))
    >>> x, y = mesh3.cellCenters
    >>> xyVar = CellVariable(mesh=mesh3, name="x y", value=x * y)
    >>> viewer3 = MatplotlibViewer(vars=numerix.sin(k * xyVar), 
    ...                            limits={'ymin': 0.1, 'ymax': 0.9}, 
    ...                            datamin=-0.9, datamax=2.0,
    ...                            title="Irregular 2D test",
    ...                            axes=ax3,
    ...                            cmap = pylab.cm.OrRd)

    >>> viewer = MultiViewer(viewers=(viewer1, viewer2, viewer3))
    >>> for kval in range(10):
    ...     k.setValue(kval)
    ...     viewer.plot()

    >>> viewer._promptForOpinion()

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """
    if type(vars) not in [type([]), type(())]:
        vars = [vars]
        
    kwlimits.update(limits)
    
    from fipy.viewers import MeshDimensionError
    from fipy.viewers.matplotlibViewer import (Matplotlib1DViewer, 
                                               Matplotlib2DGridViewer,
                                               Matplotlib2DViewer,
                                               MatplotlibVectorViewer)

    try:
        return Matplotlib1DViewer(vars=vars, title=title, axes=axes, **kwlimits)
    except MeshDimensionError:
        try:
            from fipy.viewers.matplotlibViewer.matplotlib2DGridViewer import Matplotlib2DGridViewer
            return Matplotlib2DGridViewer(vars=vars, title=title, cmap=cmap, colorbar=colorbar, axes=axes, **kwlimits)
        except MeshDimensionError:
            try:
                from fipy.viewers.matplotlibViewer.matplotlib2DViewer import Matplotlib2DViewer
                return Matplotlib2DViewer(vars=vars, title=title, cmap=cmap, colorbar=colorbar, axes=axes, **kwlimits)
            except MeshDimensionError:
                from fipy.viewers.matplotlibViewer.matplotlibVectorViewer import MatplotlibVectorViewer
                return MatplotlibVectorViewer(vars=vars, title=title, axes=axes, **kwlimits)

class _MatplotlibViewer(_Viewer):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    The `_MatplotlibViewer` is the base class for the viewers that use the
    Matplotlib_ python plotting package.
    
    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """
        
    def __init__(self, vars, title=None, figaspect=1.0, cmap=None, colorbar=None, axes=None, log=False, **kwlimits):
        """
        Create a `_MatplotlibViewer`.
        
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
        if self.__class__ is _MatplotlibViewer:
            raise NotImplementedError, "can't instantiate abstract base class"
            
        _Viewer.__init__(self, vars=vars, title=title, **kwlimits)

        import pylab

        pylab.ion()

        if axes is None:
            w, h = pylab.figaspect(figaspect)
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

        pylab.ioff()
        
        self._plot()
        pylab.draw()

        try:
            fig.canvas.flush_events()
        except NotImplementedError:
            pass
        
        pylab.ion()

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
        self._cb = matplotlib.colorbar.ColorbarBase(cbax, cmap=viewer.cmap,
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

