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

import Numeric
import Gnuplot

class GnuplotViewer:
    
    def __init__(self, array, dx = 1., dy = 1., maxVal = None, minVal = None, palette = 'color'):
        self.array = array
        self.dx = dx
        self.dy = dy
        self.maxVal = maxVal
        self.minVal = minVal
        self.palette = palette

    def plot(self, fileName = None):
        array = Numeric.array(self.array)
        (nx, ny) = array.shape
        x = Numeric.arange(ny) * self.dy
        y = Numeric.arange(nx) * self.dx
        array = Numeric.transpose(array)
        g = Gnuplot.Gnuplot()
        g('set data style lines')
##        g('set cntrparam levels auto ' + str(self.numberOfContours))
        g('set view map')
        g('set samples 101')
        g('set isosamples 5, 5')
        g('set style data pm3d')
        g('set style function pm3d')
        g('set noytics')
        g('set noztics')
        g('set noxtics')
        g('set pm3d at b')
##        g('set palette color positive')
        g('set palette ' + self.palette)
        g('set size ratio -1')
        
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

        

##        g('set zrange [ ' + str(minVal) + ' : ' + str(maxVal) + ' ]')
        g('set cbrange [ ' + str(minVal) + ' : ' + str(maxVal) + ' ]')

        if fileName is None:
            g.splot(Gnuplot.GridData(array, x, y))
        elif '.pdf' == fileName[-4:]:
            g('set terminal pdf')
            g('set output "' + fileName + '"')
            g.splot(Gnuplot.GridData(array, x, y))
        elif '.ps' == fileName[-3:]:
            g.splot(Gnuplot.GridData(array, x, y))
            g.hardcopy(fileName, enhanced=1, color=1)
