#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "matplotlib1DViewer.py"
 #                                    created: 9/14/04 {2:48:25 PM} 
 #                                last update: 5/15/06 {4:01:27 PM} { 2:45:36 PM}
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
from matplotlibViewer import MatplotlibViewer

class Matplotlib1DViewer(MatplotlibViewer):
    """
    Displays a y vs.  x plot of one or more 1D `CellVariable` objects using
    Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/


    """
    def __init__(self, vars, limits = None, title = None):
        MatplotlibViewer.__init__(self, vars=vars, limits=limits, title=title)
    
        self.lines = [pylab.plot(*datum) for datum in self._getData()]

        pylab.legend([var.getName() for var in self.vars])

        pylab.xlim(xmin = self._getLimit('xmin'),
                   xmax = self._getLimit('xmax'))

    def _getData(self):
        from fipy.tools.numerix import array
        return [[array(var.getMesh().getCellCenters()[...,0]), array(var)] for var in self.vars]
            
    def _getSuitableVars(self, vars):
        vars = [var for var in MatplotlibViewer._getSuitableVars(self, vars) if var.getMesh().getDim() == 1]
        if len(vars) > 1:
            vars = [var for var in vars if var.getMesh() is vars[0].getMesh()]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "Can only plot 1D data"
        return vars

    def _plot(self):

        ymin = self._getLimit('datamin')
        if ymin is None:
            ymin = self._getLimit('ymin')
        if ymin is None:
            ymin = min(self.vars[0])
            for var in self.vars[1:]:
                ymin = min(ymin, min(var))

        pylab.ylim(ymin=ymin)

        ymax = self._getLimit('datamax')
        if ymax is None:
            ymax = self._getLimit('ymax')
        if ymax is None:
            ymax = max(self.vars[0])
            for var in self.vars[1:]:
                ymax = max(ymax, max(var))

        pylab.ylim(ymax=ymax)

        for line, datum in zip(self.lines, self._getData()):
            line[0].set_xdata(datum[0])
            line[0].set_ydata(datum[1])
            
