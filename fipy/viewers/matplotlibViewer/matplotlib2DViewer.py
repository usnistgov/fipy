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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##
 
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from matplotlibViewer import MatplotlibViewer

class Matplotlib2DViewer(MatplotlibViewer):
    """
    Displays a contour plot of a 2D `CellVariable` object.    

    The `Matplotlib2DViewer` plots a 2D `CellVariable` using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """ 
    
    __doc__ += MatplotlibViewer._test2Dirregular(viewer="Matplotlib2DViewer")

    def __init__(self, vars, limits = None, title = None):
        """Creates a `Matplotlib2DViewer`.
        

        :Parameters:
          - `vars`: A `CellVariable` object.
          - `limits`: A dictionary with possible keys `'xmin'`, `'xmax'`, 
            `'ymin'`, `'ymax'`, `'datamin'`, `'datamax'`. Any limit set to 
            a (default) value of `None` will autoscale.
          - `title`: displayed at the top of the Viewer window

        """
        MatplotlibViewer.__init__(self, vars=vars, limits=limits, title=title, figaspect=1. / 1.3)

        self.colorbar = None
        
        self.mesh = self.vars[0].getMesh()
        
        vertexIDs = self.mesh._getOrderedCellVertexIDs()

        vertexCoords = self.mesh.getVertexCoords()

        xCoords = numerix.take(vertexCoords[0], vertexIDs)
        yCoords = numerix.take(vertexCoords[1], vertexIDs)
        
        polys = []
##         for x, y in zip(xCoords.swapaxes(0,1), yCoords.swapaxes(0,1)):
##             if hasattr(x, 'mask'):
##                 x = x.compressed()
##             if hasattr(y, 'mask'):
##                 y = y.compressed()
##             polys.append(x)
##             polys.append(y)
##             polys.append('b')

        for x, y in zip(xCoords.swapaxes(0,1), yCoords.swapaxes(0,1)):
            if hasattr(x, 'mask'):
                x = x.compressed()
            if hasattr(y, 'mask'):
                y = y.compressed()
            polys.append(zip(x,y))

        import pylab
        import matplotlib

        fig = pylab.figure(self.id)

        ax = fig.get_axes()[0]

        from matplotlib.collections import PolyCollection
        self.collection = PolyCollection(polys)
        self.collection.set_linewidth(0)
        ax.add_patch(self.collection)

        if self._getLimit('xmin') is None:
            xmin = xCoords.min()
        else:
            xmin = self._getLimit('xmin')

        if self._getLimit('xmax') is None:
            xmax = xCoords.max()
        else:
            xmax = self._getLimit('xmax')

        if self._getLimit('ymin') is None:
            ymin = yCoords.min()
        else:
            ymin = self._getLimit('ymin')

        if self._getLimit('ymax') is None:
            ymax = yCoords.max()
        else:
            ymax = self._getLimit('ymax')

        pylab.axis((xmin, xmax, ymin, ymax))

##        self.polygons = ax.fill(linewidth=0., *polys)
        
        cbax, kw = matplotlib.colorbar.make_axes(ax, orientation='vertical')
        
        # Set the colormap and norm to correspond to the data for which
        # the colorbar will be used.
        cmap = matplotlib.cm.jet
        norm = matplotlib.colors.normalize(vmin=-1, vmax=1)
        
        # ColorbarBase derives from ScalarMappable and puts a colorbar
        # in a specified axes, so it has everything needed for a
        # standalone colorbar.  There are many more kwargs, but the
        # following gives a basic continuous colorbar with ticks
        # and labels.
        self.cb = matplotlib.colorbar.ColorbarBase(cbax, cmap=cmap,
                                                   norm=norm,
                                                   orientation='vertical')
        self.cb.set_label(self.vars[0].name)
        
        self._plot()
        
    def _getSuitableVars(self, vars):
        from fipy.meshes.numMesh.mesh2D import Mesh2D
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in MatplotlibViewer._getSuitableVars(self, vars) \
          if ((isinstance(var.getMesh(), Mesh2D) and isinstance(var, CellVariable))
              and var.getRank() == 0)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "The mesh must be a Mesh2D instance"
        # this viewer can only display one variable
        return [vars[0]]
        
    def _plot(self):
##         pylab.clf()
        
##         ## Added garbage collection since matplotlib objects seem to hang
##         ## around and accumulate.
##         import gc
##         gc.collect()

        Z = self.vars[0].getValue() 
        
        zmin, zmax = self._autoscale(vars=self.vars,
                                     datamin=self._getLimit(('datamin', 'zmin')),
                                     datamax=self._getLimit(('datamax', 'zmax')))

        diff = zmax - zmin
        
        import pylab
        import matplotlib

        faceColors = []

        for value in Z:
            if diff == 0:
                rgba = pylab.cm.jet(0.5)
            else:
                rgba = pylab.cm.jet((value - zmin) / diff)

            faceColors.append(rgba)

        self.collection.set_facecolors(faceColors)
##        self.collection.set_edgecolors(faceColors)
            
##         for poly, value in zip(self.polygons, Z):
##             if diff == 0:
##                 rgba = pylab.cm.jet(0.5)
##             else:
##                 rgba = pylab.cm.jet((value - zmin) / diff)

##             poly.set_facecolor(rgba)
            
        self.cb.norm = matplotlib.colors.normalize(vmin=zmin, vmax=zmax)
        self.cb.draw_all()
        
##        pylab.xlim(xmin=self._getLimit('xmin'),
##                   xmax=self._getLimit('xmax'))

##        pylab.ylim(ymin=self._getLimit('ymin'),
##                   ymax=self._getLimit('ymax'))

    def plotMesh(self, filename = None):
        pass

if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()

        
