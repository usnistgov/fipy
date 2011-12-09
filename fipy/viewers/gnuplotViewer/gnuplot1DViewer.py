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

from fipy.variables.faceVariable import FaceVariable
from fipy.tools import numerix

from fipy.viewers.gnuplotViewer.gnuplotViewer import _GnuplotViewer
from fipy.tools.decorators import deprecateGnuplot

__all__ = ["Gnuplot1DViewer"]

@deprecateGnuplot
class Gnuplot1DViewer(_GnuplotViewer):
    """Displays a y vs. x plot of one or more 1D `CellVariable` objects.

    The `Gnuplot1DViewer` plots a 1D `CellVariable` using a front end python
    wrapper available to download (Gnuplot.py_).
    
    .. _Gnuplot.py: http://gnuplot-py.sourceforge.net/
    """
    
    __doc__ += _GnuplotViewer._test1D(viewer="Gnuplot1DViewer")

    __doc__ += """
    Different style script demos_ are available at the Gnuplot_ site.

    .. _Gnuplot: http://gnuplot.sourceforge.net/
    .. _demos: http://gnuplot.sourceforge.net/demo/

    .. note::
    
        `Gnuplot1DViewer` requires Gnuplot_ version 4.0.

    """
    def _plot(self):

        self.g('set yrange [' + self._getLimit(('datamin', 'ymin', 'zmin'))  + ':' + self._getLimit(('datamax', 'ymin', 'zmax')) + ']')
        
        tupleOfGnuplotData = ()

        import Gnuplot
        import re
        m = re.match(r"\d+.\d+", Gnuplot.__version__)
        if m is None or float(m.group(0)) < 1.8:
            raise ImportError("Gnuplot.py version 1.8 or newer is required.")
            
        for var in self.vars:
            # Python 2.6 made 'with' a keyward (deprecation warnings have been issued since 2.5)
            # this was addressed in Gnuplot.py in r299, in 2007
            
            if var._variableClass is not FaceVariable:
                X = var.mesh.cellCenters[0]
            else:
                X = var.mesh.faceCenters[0]

            tupleOfGnuplotData += (Gnuplot.Data(numerix.array(X), 
                                                numerix.array(var),
                                                title=var.name,
                                                with_='lines'),)
                              
        apply(self.g.plot, tupleOfGnuplotData)
    
if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
