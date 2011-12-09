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

from fipy.meshes import Grid2D
from fipy.variables.faceVariable import FaceVariable
from fipy.tools import numerix

from fipy.viewers.gnuplotViewer.gnuplotViewer import _GnuplotViewer
from fipy.tools.decorators import deprecateGnuplot

__all__ = ["Gnuplot2DViewer"]

@deprecateGnuplot
class Gnuplot2DViewer(_GnuplotViewer):
    """
    Displays a contour plot of a 2D `CellVariable` object.    
       
    The `Gnuplot2DViewer` plots a 2D `CellVariable` using a front end
    python wrapper available to download (Gnuplot.py_).

    .. _Gnuplot.py: http://gnuplot-py.sourceforge.net/
    """
    
    __doc__ += _GnuplotViewer._test2D(viewer="Gnuplot2DViewer")
    __doc__ += _GnuplotViewer._test2Dirregular(viewer="Gnuplot2DViewer")
    
    __doc__ += """
    Different style script demos_ are available at the Gnuplot_ site.

    .. _Gnuplot: http://gnuplot.sourceforge.net/
    .. _demos: http://gnuplot.sourceforge.net/demo/

    .. note::
    
        `Gnuplot2DViewer` requires Gnuplot_ version 4.0.

    """
    def __init__(self, vars, title = None, limits={}, **kwlimits):
        """Creates a `Gnuplot2DViewer`.

        :Parameters:
          vars
            a `CellVariable` object.
          title
            displayed at the top of the `Viewer` window
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to 
            a (default) value of `None` will autoscale.
        """
        kwlimits.update(limits)
        _GnuplotViewer.__init__(self, vars=vars, title=title, **kwlimits)
        
        if len(self.vars) != 1:
            raise IndexError, "A 2D Gnuplot viewer can only display one Variable"
            
    def _plot(self):

        self.g('set cbrange [' + self._getLimit(('datamin', 'zmin'))  + ':' + self._getLimit(('datamax', 'zmax')) + ']')
        self.g('set view map')
        self.g('set style data pm3d')
        self.g('set pm3d at st solid')
        mesh = self.vars[0].mesh

        if self.vars[0]._variableClass is FaceVariable:
            x, y = mesh.faceCenters
            if isinstance(mesh, Grid2D.__class__):
                nx, ny = mesh.shape
            else:
                N = int(numerix.sqrt(mesh.numberOfCells))
                nx, ny = N, N
        else:
            x, y = mesh.cellCenters
            N = int(numerix.sqrt(mesh.numberOfFaces))
            nx, ny = N, N

        self.g('set dgrid3d %i, %i, 2' % (ny, nx))

        import Gnuplot
        data = Gnuplot.Data(numerix.array(x), numerix.array(y),
                            self.vars[0].value)

        self.g.splot(data)

if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
