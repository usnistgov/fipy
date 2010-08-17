#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "matplotlib2DViewer.py"
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

from matplotlibViewer import _MatplotlibViewer, _ColorBar

class Matplotlib2DGridViewer(_MatplotlibViewer):
    """
    Displays an image plot of a 2D `CellVariable` object using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """
    
    __doc__ += _MatplotlibViewer._test2D(viewer="Matplotlib2DGridViewer")

    def __init__(self, vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, **kwlimits):
        """
        Creates a `Matplotlib2DGridViewer`.
        
        :Parameters:
          vars
            A `CellVariable` object.
          title
            displayed at the top of the `Viewer` window
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          cmap
            The colormap. Defaults to `matplotlib.cm.jet`
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to 
            a (default) value of `None` will autoscale.
          colorbar
            plot a colorbar in specified orientation if not `None`
          axes
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        """
        kwlimits.update(limits)
        _MatplotlibViewer.__init__(self, vars=vars, title=title, 
                                   cmap=cmap, colorbar=colorbar, axes=axes, 
                                   **kwlimits)

        self.image = self.axes.imshow(self._getData(),
                                      extent=(self._getLimit('xmin'), self._getLimit('xmax'), 
                                              self._getLimit('ymin'), self._getLimit('ymax')),
                                      vmin=self._getLimit(key=('datamin', 'zmin')),
                                      vmax=self._getLimit(key=('datamax', 'zmax')),
                                      cmap=self.cmap)
                   
        if title is None:                          
            self.axes.set_title(self.vars[0].getName())

    def _getLimit(self, key, default=None):
        limit = _MatplotlibViewer._getLimit(self, key, default=default)
        if limit is None:
            X, Y = self.vars[0].getMesh().getFaceCenters()
            if 'xmin' in key:
                limit = float(min(X))
            elif 'ymin' in key:
                limit = float(min(Y))
            if 'xmax' in key:
                limit = float(max(X))
            elif 'ymax' in key:
                limit = float(max(Y))
        return limit
        
    def _getSuitableVars(self, vars):
        from fipy.meshes.numMesh.uniformGrid2D import UniformGrid2D
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in _MatplotlibViewer._getSuitableVars(self, vars) \
          if (isinstance(var.getMesh(), UniformGrid2D) and isinstance(var, CellVariable)
              and var.getRank() == 0)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "Matplotlib2DGridViewer can only display a rank-0 CellVariable with a UniformGrid2D mesh"
        # this viewer can only display one variable
        return [vars[0]]
        
    def _getData(self):
        from fipy.tools.numerix import array, reshape
        return reshape(array(self.vars[0]), self.vars[0].getMesh().getShape()[::-1])[::-1]

    def _plot(self):
        datamin = self._getLimit(('datamin', 'zmin')) 
        datamax = self._getLimit(('datamax', 'zmax')) 
        if datamin is None or datamax is None:
            datamin, datamax = self._autoscale(vars=self.vars, 
                                               datamin=datamin, 
                                               datamax=datamax)

            self.image.set_clim(vmax=datamax, vmin=datamin)

        self.image.set_data(self._getData())
        
        if self.colorbar is not None:
            self.colorbar.plot(vmin=datamin, vmax=datamax)

def _test():
    from fipy.viewers.viewer import _test2D
    _test2D(Matplotlib2DGridViewer)

if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
