#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "matplotlibVectorViewer.py"
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
from matplotlibViewer import _MatplotlibViewer
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable

class MatplotlibVectorViewer(_MatplotlibViewer):
    """Displays a vector plot of a 2D rank-1 `CellVariable` or
    `FaceVariable` object using Matplotlib_

    .. _Matplotlib: http://matplotlib.sourceforge.net/

    """
    
    __doc__ += _MatplotlibViewer._test2Dvector(viewer="MatplotlibVectorViewer")
    __doc__ += """
    
            >>> for sparsity in arange(5000, 0, -500):
            ...     viewer.quiver(sparsity=sparsity)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        
    """
    __doc__ += _MatplotlibViewer._test2DvectorIrregular(viewer="MatplotlibVectorViewer")

    def __init__(self, vars, title=None, scale=None, sparsity=None, limits={}, axes=None, **kwlimits):
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
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to 
            a (default) value of `None` will autoscale.
          axes
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        """
        kwlimits.update(limits)
        _MatplotlibViewer.__init__(self, vars=vars, title=title, axs=axes, **kwlimits)

        self.quiver(sparsity=sparsity, scale=scale)
        
        self._plot()
        
    def quiver(self, sparsity=None, scale=None):
        var = self.vars[0]
        mesh = var.getMesh()

        if isinstance(var, FaceVariable):
            N = mesh._getNumberOfFaces() 
            V = mesh._getFaceAreas()
            X, Y = mesh.getFaceCenters()
        elif isinstance(var, CellVariable):
            N = mesh.getNumberOfCells() 
            V = mesh.getCellVolumes()
            X, Y = mesh.getCellCenters()

        if sparsity is not None and N > sparsity:
            self.indices = numerix.random.rand(N) * V
            self.indices = self.indices.argsort()[-sparsity:]
        else:
            self.indices = numerix.arange(N)

        X = numerix.take(X, self.indices)
        Y = numerix.take(Y, self.indices)
        
        U = V = numerix.ones(X.shape)
        
        import pylab
        
        pylab.ion()
        self.axes.cla()
        self._quiver = self.axes.quiver(X, Y, U, V, scale=scale)
        pylab.ioff()

    def _getSuitableVars(self, vars):
        from fipy.meshes.numMesh.mesh2D import Mesh2D

        vars = [var for var in _MatplotlibViewer._getSuitableVars(self, vars) \
                if (isinstance(var.getMesh(), Mesh2D) \
                    and (isinstance(var, FaceVariable) \
                         or isinstance(var, CellVariable)) and var.getRank() == 1)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "The mesh must be a Mesh2D instance"
        # this viewer can only display one variable
        return [vars[0]]
                
    def _plot(self):

        var = self.vars[0]
        mesh = var.getMesh()

        U, V = var.getNumericValue()

        U = numerix.take(U, self.indices)
        V = numerix.take(V, self.indices)

        self._quiver.set_UVC(U, V)
        
        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))
        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
