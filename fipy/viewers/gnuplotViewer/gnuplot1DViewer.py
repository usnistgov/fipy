#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gnuplot1DViewer.py"
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

"""

The `Gnuplot1DViewer` plots a 1D Numeric array using a front end python
wrapper available to download (Gnuplot.py_).

.. _Gnuplot.py: http://gnuplot-py.sourceforge.net/

Different style script demos_ are available at the Gnuplot_ site.

.. _Gnuplot: http://gnuplot.sourceforge.net/
.. _demos: http://gnuplot.sourceforge.net/demo/

.. note::
    
   `GnuplotViewer` requires Gnuplot_ version 4.0.

"""
__docformat__ = 'restructuredtext'

import Numeric
import Gnuplot

from gnuplotViewer import GnuplotViewer

class Gnuplot1DViewer(GnuplotViewer):
    
    def _plot(self):
        
        mesh = self.vars[0].getMesh()
        x = mesh.getCellCenters()[:,0]
        NCells = mesh.getNumberOfCells()
##        listOfArrays = []
        listOfGnuplotData = ()

        for var in self.vars:
##            arr = Numeric.zeros((NCells, 2), 'd')
##            arr[:,0] = x
##            arr[:,1] = var[:]
##            listOfArrays += [arr]
            listOfGnuplotData += (Gnuplot.Data(mesh.getCellCenters()[:,0],
                                               var[:],
                                               title=var.getName(),
                                               with='lines'),)

            
        self.g('set data style linespoints')
        apply(self.g.plot, listOfGnuplotData)
    
