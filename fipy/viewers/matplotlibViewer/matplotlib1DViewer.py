#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "matplotlib1DViewer.py"
 #                                    created: 9/14/04 {2:48:25 PM} 
 #                                last update: 9/2/05 {10:48:31 AM} { 2:45:36 PM}
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
    
    def _plot(self):

        data = ()
        names = ()
        for var in self.vars:
            data += (Numeric.array(var.getMesh().getCellCenters()[:,0]), Numeric.array(var))
            names += (var.getName(),)

        apply(pylab.plot, data)
        pylab.legend(names)
        
        ymin = self._getLimit('datamin')
        if ymin is None:
            ymin = self._getLimit('ymin')
        pylab.ylim(ymin = ymin)

        ymax = self._getLimit('datamax')
        if ymax is None:
            ymax = self._getLimit('ymax')
        pylab.ylim(ymax = ymax)
        
