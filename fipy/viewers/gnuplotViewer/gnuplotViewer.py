#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gnuplotViewer.py"
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

The `GnuplotViewer` is the base class for `Gnuplot1DViewer` and
`Gnuplot2DViewer` It using a front end python wrapper available to
download (Gnuplot.py_).

.. _Gnuplot.py: http://gnuplot-py.sourceforge.net/

If one would like more specific styles for a plot, it is probably best
to create an inherited class and change the `gnuplotCommands` method.

Different style script demos_ are available at the Gnuplot_ site.

.. _Gnuplot: http://gnuplot.sourceforge.net/
.. _demos: http://gnuplot.sourceforge.net/demo/

.. note::
    
   `GnuplotViewer` requires Gnuplot_ version 4.0.

"""
__docformat__ = 'restructuredtext'

import Gnuplot
from fipy.viewers.viewer import Viewer

class GnuplotViewer(Viewer):
    
    def __init__(self, vars, limits = None, title = None):
        """

        :Parameters:

          - `vars`: a `Variable` or tuple of `Variable` objects to plot
          - `limits`: a dictionary with possible keys `xmin`, `xmax`,
                      `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.
                      A 1D Viewer will only use `xmin` and `xmax`, a 2D viewer
                      will also use `ymin` and `ymax`, and so on.
                      All viewers will use `datamin` and `datamax`.
                      Any limit set to a (default) value of `None` will autoscale.
          - `title`: displayed at the top of the Viewer window

        """
        Viewer.__init__(self, vars = vars, limits = limits, title = title)
        self.g = Gnuplot.Gnuplot()
        self.g('set title "' + self.title + '"')

    def getLimit(self, key):
        limit = Viewer.getLimit(self, key)
        if limit is None:
            return ''
        else:
            return str(limit)

    def plot(self, fileName = None):

        self.g('set xrange [' + self.getLimit('xmin')  + ':' + self.getLimit('xmax') + ']')
        self.g('set yrange [' + self.getLimit('ymin')  + ':' + self.getLimit('ymax') + ']')
        self.g('set zrange [' + self.getLimit('zmin')  + ':' + self.getLimit('zmax') + ']')
        self.g('set cbrange [' + self.getLimit('datamin')  + ':' + self.getLimit('datamax') + ']')

        self._plot()
        if fileName is not None:
            self.g.hardcopy(fileName)

