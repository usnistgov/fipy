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
    
    def __init__(self, array, numberOfContours = 10, dx = 1., dy = 1.):
        self.array = array
        self.numberOfContours = numberOfContours
        self.dx = dx
        self.dy = dy

    def plot(self, fileName = None):
        array = Numeric.array(self.array)
        (nx, ny) = array.shape
        x = Numeric.arange(ny) * self.dy
        y = Numeric.arange(nx) * self.dx
        array = Numeric.transpose(array)
        g = Gnuplot.Gnuplot(debug=1)
        g('set data style lines')
##        g('set contour base')
##        g('set nosurface')
##        g('set view map')
        g('set cntrparam levels auto ' + str(self.numberOfContours))
##        g('set view map')
##        g('unset surface')
##        g('set size ratio -1')
##        g('unset xtics')
##        g('unset ytics')
##        g('set palette model RGB')
##        g('set palette defined ( 0 0.05 0.05 0.2, 0.1 0 0 1, 0.25 0.7 0.85 0.9, 0.4 0 0.75 0, 0.5 1 1 0, 0.7 1 0 0, 0.9 0.6 0.6 0.6, 1 0.95 0.95 0.95 )')
        g('set view map')
        g('set samples 101')
        g('set isosamples 5, 5')
        g('set style data pm3d')
        g('set style function pm3d')
        g('set noytics')
        g('set noztics')
        g('set noxtics')
##        g('set cbrange [ -10.0000 : 10.0000 ] noreverse nowriteback')
        g('set pm3d at b')
        g('set palette defined ( 0 0.05 0.05 0.2, 0.1 0 0 1, 0.25 0.7 0.85 0.9, 0.4 0 0.75 0, 0.5 1 1 0, 0.7 1 0 0, 0.9 0.6 0.6 0.6, 1 0.95 0.95 0.95 )')
        g('set pm3d at b')
        g.splot(Gnuplot.GridData(array, x, y))

        if fileName is not None:

            if '.pdf' == fileName[-4:]:
                g('set terminal pdf')
                g('set output "' + fileName + '"')
                g.splot(Gnuplot.GridData(array, x, y))

            if '.ps' == fileName[-3:]: 
                g.hardcopy(fileName, enhanced=1, color=1)
