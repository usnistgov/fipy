#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gnuplotViewer.py"
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

The `GnuplotViewer` plots a 2D Numeric array using a front end python
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

class GnuplotViewer:
    
    def __init__(self, array, dx = 1., dy = 1., maxVal = None, minVal = None, palette = 'color', legend = True, title = ''):
        """

        Argument list:

        `array` - A 2D Numeric array

        `dx` - The unit size is the x-direction.

        `dy` - The unit size in the y-direction.

        `maxVal` - The maximum value to appear on the legend

        `maxVal` - The minimum value to appear on the legend.

        `palette` - The color scheme. Other choices might be `color
        negative` or `gray`. This is completely configurable.

        `legend` - Whether to display the legend.
        
        """
        
        self.array = array
        self.dx = dx
        self.dy = dy
        self.maxVal = maxVal
        self.minVal = minVal
        self.palette = palette
        self.legend = legend
        self.title = title
        
    def gnuplotCommands(self, g):
        g('set data style lines')
        g('set view map')
        g('set samples 101')
        g('set isosamples 5, 5')
        g('set style data pm3d')
        g('set style function pm3d')
        g('set noytics')
        g('set noztics')
        g('set noxtics')
        g('set pm3d at b')
        g('set palette ' + self.palette)
        g('set size ratio -1')
        g('set title "' + self.title + '"')
##        g('set size square')
        return g

    def plot(self, fileName = None):
        """

        Argument list:

        `fileName` - If no file name is passed a plot appears on the
        screen.  If `*.ps` or `*.pdf` is passed, a postscript or PDF
        file is produced and nothing appears on the screen.

        """
        
        ## prepare the array
        array =  Numeric.array(self.array)
        (nx, ny) = array.shape
        x = Numeric.arange(ny) * self.dy
        y = Numeric.arange(nx) * self.dx
        array = Numeric.transpose(array)

        ## issue the gnuplot commands

        g = Gnuplot.Gnuplot()
        g = self.gnuplotCommands(g)

        ## obtain the max and min values
        
        if self.maxVal is not None:
            maxVal = self.maxVal
        else:
            argmax = Numeric.argmax(array.flat)
            maxVal = array.flat[argmax]

        if self.minVal is not None:
            minVal = self.minVal
        else:
            argmin = Numeric.argmin(array.flat)
            minVal = array.flat[argmin]

        ## write to the screen or a file.

##        g('set zrange [ ' + str(minVal) + ' : ' + str(maxVal) + ' ]')
        g('set cbrange [ ' + str(minVal) + ' : ' + str(maxVal) + ' ]')

        if not self.legend:
            g('unset colorbox')
            g('unset border')
                
        if fileName is None:
            g.splot(Gnuplot.GridData(array, x, y))
        elif '.pdf' == fileName[-4:]:
            g('set terminal pdf')
            g('set output "' + fileName + '"')
            g.splot(Gnuplot.GridData(array, x, y))
        elif '.ps' == fileName[-3:]:
            g.splot(Gnuplot.GridData(array, x, y))
            g.hardcopy(fileName, enhanced=1, color=1)
