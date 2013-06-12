#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "matplotlibStreamViewer.py"
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

from fipy.tools import numerix
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable

from fipy.viewers.matplotlibViewer.matplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["MatplotlibStreamViewer"]

class MatplotlibStreamViewer(AbstractMatplotlib2DViewer):
    """Displays a stream plot of a 2D rank-1 `CellVariable` or
    `FaceVariable` object using Matplotlib_

    One issue is that this `Viewer` relies on `scipy.interpolate.griddata`,
    which interpolates on the convex hull of the data. The results is that
    streams are plotted across any concavities in the mesh.
    
    .. _Matplotlib: http://matplotlib.sourceforge.net/

    """
    
#     __doc__ += AbstractMatplotlib2DViewer._test2Dvector(viewer="MatplotlibStreamViewer")
#     __doc__XXX += """
#     
#            >>> for sparsity in numerix.arange(5000, 0, -500):
#            ...     viewer.stream(sparsity=sparsity)
#            ...     viewer.plot()
#            >>> viewer._promptForOpinion()
#         
#     """
    __doc__ += AbstractMatplotlib2DViewer._test2DvectorIrregular(viewer="MatplotlibStreamViewer")

    def __init__(self, vars, title=None, scale=None, sparsity=None, log=False, limits={}, axes=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DViewer`.

        :Parameters:
          vars
            a rank-1 `CellVariable` or `FaceVariable` object.
          title
            displayed at the top of the `Viewer` window
          scale
            if not `None`, scale all arrow lengths by this value
          sparsity
            if not `None`, then this number of arrows will be
            randomly chosen (weighted by the cell volume or face area)
          log
            if `True`, arrow length goes at the base-10 logarithm of the magnitude
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to 
            a (default) value of `None` will autoscale.
          axes
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
          figaspect
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is 'auto', the aspect ratio will be determined from
            the Variable's mesh.
        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title, axes=axes, figaspect=figaspect, **kwlimits)

        self.log = log
        
        self._plot()
        
    def _getSuitableVars(self, vars):
        from fipy.meshes.mesh2D import Mesh2D
        from fipy.meshes.uniformGrid2D import UniformGrid2D

        vars = [var for var in AbstractMatplotlib2DViewer._getSuitableVars(self, vars) \
                if ((isinstance(var.mesh, Mesh2D) 
                     or isinstance(var.mesh, UniformGrid2D))\
                    and (isinstance(var, FaceVariable) \
                         or isinstance(var, CellVariable)) and var.rank == 1)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "The mesh must be a Mesh2D instance"
        # this viewer can only display one variable
        return [vars[0]]
                
    def _plot(self):
        from scipy.interpolate import griddata

        var = self.vars[0]
        mesh = var.mesh

        vertexIDs = mesh._orderedCellVertexIDs

        vertexCoords = mesh.getVertexCoords()

        X = numerix.take(vertexCoords[0], vertexIDs)
        Y = numerix.take(vertexCoords[1], vertexIDs)

        xmin = X.min()
        xmax = X.max()
        ymin = Y.min()
        ymax = Y.max()
        
        N = 100
        X = numerix.linspace(xmin, xmax, N)
        Y = numerix.linspace(ymin, ymax, N)
                                
        grid_x, grid_y = numerix.mgrid[xmin:xmax:N*1j, ymin:ymax:N*1j]
        
        U = griddata(mesh.cellCenters.value.T,
                     var.value[0], (grid_x, grid_y), method='cubic')
        V = griddata(mesh.cellCenters.value.T, 
                     var.value[1], (grid_x, grid_y), method='cubic')
                     
        U = U.T
        V = V.T

        ang = numerix.arctan2(V, U)
        mag = numerix.sqrt(U**2 + V**2)
        
        datamin, datamax = self._autoscale(vars=(mag,),
                                           datamin=self._getLimit('datamin'),
                                           datamax=self._getLimit('datamax'))
        
        mag = numerix.where(mag > datamax, numerix.nan, mag)
        mag = numerix.where(mag < datamin, numerix.nan, mag)
        
        if self.log:
            mag = numerix.log10(mag)
            
        U = mag * numerix.cos(ang)
        V = mag * numerix.sin(ang)

        self.axes.cla()
        self.axes.streamplot(X, Y, U, V)
        
        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))
        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
