#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "matplotlibViewer.py"
 #                                    created: 9/14/04 {2:48:25 PM} 
 #                                last update: 7/6/05 {4:44:23 PM} { 2:45:36 PM}
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

from fipy.viewers.viewer import Viewer

class MatplotlibViewer(Viewer):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    The `MatplotlibViewer` is the base class for the viewers that use the
    Matplotlib_ python plotting package.

    .. _Matplotlib: http://matplotlib.sourceforge.net/

    """
        
    def __init__(self, vars, limits = None, title = None):
        """
        Create a `MatplotlibViewer`.
        
        :Parameters:

          - `vars`: a `CellVariable` or tuple of `CellVariable` objects to plot
          - `limits`: a dictionary with possible keys `xmin`, `xmax`,
            `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.  A 1D
            Viewer will only use `xmin` and `xmax`, a 2D viewer will also
            use `ymin` and `ymax`, and so on.  All viewers will use
            `datamin` and `datamax`.  Any limit set to a (default) value of
            `None` will autoscale.
          - `title`: displayed at the top of the Viewer window

        """
        Viewer.__init__(self, vars = vars, limits = limits, title = title)

        pylab.ion()

        fig = pylab.figure()
        self.id = fig.number
                    
    def plot(self, fileName = None):
        """
        Plot the `CellVariable` as a contour plot.

        :Parameters:
          - `fileName`: The name of the file for hard copies.
          
        """

        pylab.figure(self.id)
        pylab.clf()
        pylab.title(self.title)
        
        self._plot()

        pylab.xlim(xmin = self._getLimit('xmin'))
        pylab.xlim(xmax = self._getLimit('xmax'))
        
        if fileName is not None:
            pylab.savefig(fileName)

