#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gnuplot1DViewer.py"
 #                                    created: 9/14/04 {2:48:25 PM} 
 #                                last update: 9/14/04 {10:30:12 PM} { 2:45:36 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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

The `GnuplotViewer` plots a 1D Numeric array using a front end python
wrapper available to download (Gnuplot.py_).

.. _Gnuplot.py: http://gnuplot-py.sourceforge.net/

If one would like more specific styles for a plot, it is probably best
to create an inherited class and change the `gnuplotCommands` method.

Different style script demos_ are available at the Gnuplot_ site.

.. _Gnuplot: http://gnuplot.sourceforge.net/
.. _demos: http://gnuplot.sourceforge.net/demo/

`GnuplotViewer` requires Gnuplot_ version 4.0.

"""
__docformat__ = 'restructuredtext'

import Numeric
import Gnuplot

class Gnuplot1DViewer:
    
    def __init__(self, vars = (), title = '', varTitles = (), xlabel = '', ylabel = '', keycoord = (1,1)):
        """

        Argument list:

        `var` - The variable to be plotted

        """
        
        self.vars = vars
        self.varTitles = varTitles
        self.plotter = Gnuplot.Gnuplot()
        self.plotter('set size .5')
        self.plotter('set key ' + str(keycoord[0]) + ',' + str(keycoord[1]))
        self.plotter('show key')
        self.plotter.xlabel(xlabel)
        self.plotter.ylabel(ylabel)
        self.plotter.title(title)
        
    def plot(self, fileName = ''):

        data = ()
        i = 0
        for var in self.vars:

            data += (Gnuplot.Data(var.getMesh().getCellCenters()[:,0],
                                  var,
                                  title = self.varTitles[i],
                                  with = 'lines lw 6 lt ' + str(i + 1)),)
            i += 1

        self.plotter.plot(data[0], data[1], data[2], data[3])

        if '.ps' == fileName[-3:]:
            self.plotter('set xtics 0.2')
            self.plotter('set ytics .5')
            self.plotter('set xtics nomirror')
            self.plotter('set ytics nomirror')
            self.plotter('set terminal postscript landscape enhanced color dashed "arial" 15')
            self.plotter('set output "' + fileName + '"')
            self.plotter.plot(data[0], data[1], data[2], data[3])
##            self.plotter.hardcopy(fileName, enhanced = 1, color = 0, fontsize = 21)

        elif '.pdf' == fileName[-4:]:
            self.plotter('set terminal pdf fsize 10')
            self.plotter('set output "' + fileName + '"')
            self.plotter.plot(data[0], data[1], data[2], data[3])
