#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gnuplotViewer.py"
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

from fipy.viewers.viewer import Viewer

class GnuplotViewer(Viewer):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    The `GnuplotViewer` is the base class for `Gnuplot1DViewer` and
    `Gnuplot2DViewer` It uses a front end python wrapper available to
    download (Gnuplot.py_).

    .. _Gnuplot.py: http://gnuplot-py.sourceforge.net/

    Different style script demos_ are available at the Gnuplot_ site.

    .. _Gnuplot: http://gnuplot.sourceforge.net/
    .. _demos: http://gnuplot.sourceforge.net/demo/

    .. note::
    
        `GnuplotViewer` requires Gnuplot_ version 4.0.

   """    
    def __init__(self, vars, limits = None, title = None):
        """
        The `GnuplotViewer` should not be called directly only `Gnuplot1DViewer`
        and `Gnuplot2DViewer` should be called.
        
        :Parameters:

          - `vars`: a `CellVariable` or tuple of `CellVariable` objects to plot
          - `limits`: a dictionary with possible keys `xmin`, `xmax`,
            `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.  A 1D
            Viewer will only use `xmin` and `xmax`, a 2D viewer will also
            use `ymin` and `ymax`, and so on.  All viewers will use
            `datamin` and `datamax`.  Any limit set to a (default) value of
            `None` will autoscale.
          - `title`: displayed at the top of the Viewer window

        """
        if self.__class__ is GnuplotViewer:
            raise NotImplementedError, "can't instantiate abstract base class"
            
        Viewer.__init__(self, vars = vars, limits = limits, title = title)
        import Gnuplot
        self.g = Gnuplot.Gnuplot()
        self.g('set title "' + self.title + '"')

    def _getLimit(self, key):
        limit = Viewer._getLimit(self, key)
        if limit is None:
            return ''
        else:
            return str(limit)

    def plot(self, filename = None):
        """
        Plot the `CellVariable` as a contour plot.

        :Parameters:
          - `filename`: The name of the file for hard copies.
          
        """
        pairs = (('x', 'x'), ('y', 'y'), ('z', 'z'), ('cb', 'data'))
        
        for pair  in pairs:
            self.g('set ' + pair[0] + 'range [' + self._getLimit(pair[1] + 'min')  + ':' + self._getLimit(pair[1] + 'max') + ']')

        self._plot()
        if filename is not None:
            self.g.hardcopy(filename)

