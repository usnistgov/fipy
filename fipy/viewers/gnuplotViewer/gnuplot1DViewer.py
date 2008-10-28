#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gnuplot1DViewer.py"
 #
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
 # ###################################################################
 ##
 
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from gnuplotViewer import GnuplotViewer

class Gnuplot1DViewer(GnuplotViewer):
    """Displays a y vs. x plot of one or more 1D `CellVariable` objects.

    The `Gnuplot1DViewer` plots a 1D `CellVariable` using a front end python
    wrapper available to download (Gnuplot.py_).
    
    .. _Gnuplot.py: http://gnuplot-py.sourceforge.net/
    """
    
    __doc__ += GnuplotViewer._test1D(viewer="Gnuplot1DViewer")

    __doc__ += """
    Different style script demos_ are available at the Gnuplot_ site.

    .. _Gnuplot: http://gnuplot.sourceforge.net/
    .. _demos: http://gnuplot.sourceforge.net/demo/

    .. note::
    
        `GnuplotViewer` requires Gnuplot_ version 4.0.

    """
    def _plot(self):

        self.g('set yrange [' + self._getLimit(('datamin', 'ymin', 'zmin'))  + ':' + self._getLimit(('datamax', 'ymin', 'zmax')) + ']')
        
        tupleOfGnuplotData = ()

        import Gnuplot
        for var in self.vars:
            tupleOfGnuplotData += (Gnuplot.Data(numerix.array(var.getMesh().getCellCenters()[0]),
                                                numerix.array(var.getValue()),
                                                title=var.getName(),
                                                with='lines'),)

        apply(self.g.plot, tupleOfGnuplotData)
    
if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
