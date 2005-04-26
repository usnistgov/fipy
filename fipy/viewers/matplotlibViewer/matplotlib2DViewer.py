#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "matplotlib2DViewer.py"
 #                                    created: 9/14/04 {2:48:25 PM} 
 #                                last update: 11/16/04 {10:15:25 AM} { 2:45:36 PM}
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

import pylab
import Numeric
from matplotlibViewer import _MatplotlibViewer

class Matplotlib2DViewer(_MatplotlibViewer):
    """
    Displays a contour plot of a 2D `CellVariable` object.    
    Usage

    ::
    
       viewer = Matplotlib2DViewer(var)
       viewer.plot()

    The `Matplotlib2DViewer` plots a 2D `CellVariable` using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/


    """


    def __init__(self, vars, limits = None, title = None):
        """
        Creates a `Matplotlib2DViewer`.
        
        :Parameters:
          - `vars`: A `CellVariable` object.
          - `limits`: A dictionary with possible keys `'xmin'`, `'xmax'`, 
            `'ymin'`, `'ymax'`, `'zmin'`, `'zmax'`, `'datamin'`, `'datamax'`.  A 1D
            Viewer will only use `'xmin'` and `'xmax'`, a 2D viewer will also
            use `'ymin'` and `'ymax'`, and so on.  All viewers will use
            `'datamin'` and `'datamax'`.  Any limit set to a (default) value of
            `None` will autoscale.
          - `title`: displayed at the top of the Viewer window

        """
        _MatplotlibViewer.__init__(self, vars = vars, limits = limits, title = title)
        
        if len(self.vars) != 1:
            raise IndexError, "A 2D Matplotlib viewer can only display one Variable"

        from fipy.meshes.grid2D import Grid2D
        if not  isinstance(self.vars[0].getMesh(), Grid2D):
            raise 'The mesh must be a Grid2D instance for the Matpoltlib2dViewer'

        self.colorbar = False
        
    def _plot(self):

        mesh = self.vars[0].getMesh()
        shape = mesh.getShape()
        shape = (shape[1], shape[0])
        X = Numeric.reshape(mesh.getCellCenters()[:,0], shape)
        Y = Numeric.reshape(mesh.getCellCenters()[:,1], shape)
        Z = Numeric.reshape(self.vars[0][:], shape)
        
        minz = min(self.vars[0])
        for limit in ('zmin', 'datamin'):
            value = self._getLimit(limit)
            if value is not None:
                minz = max(min(self.vars[0]), value)

        maxz = max(self.vars[0])
        for limit in ('zmax', 'datamax'):
            value = self._getLimit(limit)
            if value is not None:
                maxz = min(max(self.vars[0]), value)

        numberOfContours = 10
        smallNumber = 1e-8
        diff = maxz - minz
        
        if diff < smallNumber:
            diff = smallNumber
            
        V = Numeric.arange(numberOfContours + 1) * diff / numberOfContours + minz
            
        pylab.hot()
        pylab.contour(X, Y, Numeric.reshape(self.vars[0][:], shape), V)
        ##pylab.contourf(X, Y, Numeric.reshape(self.vars[0][:], shape), V)

        ##if self.colorbar is False:
        pylab.colorbar()
        ##    self.colorbar = True
            
        pylab.ylim(ymin = self._getLimit('ymin'))
        pylab.ylim(ymax = self._getLimit('ymax'))


        
