#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mayaviDistanceViewer.py"
 #                                    created: 7/29/04 {10:39:23 AM} 
 #                                last update: 1/12/06 {7:42:58 PM}
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.viewers.viewer import Viewer

class MayaviDistanceViewer(Viewer):
    
    """
    The `MayaviDistanceViewer` creates a viewer with the Mayavi_ python
    plotting package that displays a `DistanceVariable`.

    .. _Mayavi: http://mayavi.sourceforge.net/

    """
        
    def __init__(self, distanceVar, surfactantVar, levelSetValue = 0., limits = None, title = None):
        """
        Create a `MayaviDistanceViewer`.
        
        :Parameters:

          - `distanceVar`: a `DistanceVariable` object.
          - `levelSetValue`: the value of the contour to be displayed
          - `limits`: a dictionary with possible keys `xmin`, `xmax`,
            `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.  A 1D
            Viewer will only use `xmin` and `xmax`, a 2D viewer will also
            use `ymin` and `ymax`, and so on.  All viewers will use
            `datamin` and `datamax`.  Any limit set to a (default) value of
            `None` will autoscale.
          - `title`: displayed at the top of the Viewer window
        """

        Viewer.__init__(self, vars = [], limits = limits, title = title)
        import mayavi
        self._viewer = mayavi.mayavi()
        self.distanceVar = distanceVar
        self.surfactantVar = surfactantVar
        if distanceVar.getMesh().getDim() != 2:
            raise 'The MayaviIsoViewer only works for 2D meshes.'

    def _getStructure(self):

        IDs = numerix.nonzero(self.distanceVar._getCellInterfaceFlag())
        coordinates = numerix.take(numerix.array(self.distanceVar.getMesh().getCellCenters()), IDs)
        coordinates -= numerix.take(self.distanceVar.getGrad() * self.distanceVar / self.distanceVar.getGrad().getMag(), IDs)

        from lines import _getOrderedLines
        lines = _getOrderedLines(range(len(IDs)), coordinates)
        argDict = {'points': coordinates, 'poly_line': lines}
##        for line in lines:
##            argDict['poly_line'] = line

        data = numerix.take(self.surfactantVar, IDs)

        import pyvtk
        return (pyvtk.UnstructuredGrid(points = coordinates,
                                       poly_line = lines),
                pyvtk.PointData(pyvtk.Scalars(data)))
        
    def plot(self, filename = None):

        import pyvtk
        structure, data = self._getStructure()
        data = pyvtk.VtkData(structure, data)

        import tempfile
        (f, tempFileName) = tempfile.mkstemp('.vtk')
        data.tofile(tempFileName)
        self._viewer.open_vtk(tempFileName, config=0)

        import os
        os.close(f)
        os.remove(tempFileName)
        self._viewer.load_module('SurfaceMap', 0)
        rw = self._viewer.get_render_window()
        rw.z_plus_view()

        ## display legend
        dvm = self._viewer.get_current_dvm()
        mm = dvm.get_current_module_mgr()
        slh = mm.get_scalar_lut_handler()
        slh.legend_on.set(1)
        slh.legend_on_off()
        
        ## display legen with correct range
        slh.range_on_var.set(1)
        slh.v_range_on_var.set(1)
        
        xmax = self._getLimit('datamax')
        if xmax is None:
            xmax = numerix.max(self.surfactantVar)
            
        xmin = self._getLimit('datamin')
        if xmin is None:
            xmin = numerix.min(self.surafactantVar)
            
        slh.range_var.set((xmin, xmax))
        slh.set_range_var()
        
        slh.v_range_var.set((numerix.min(var), numerix.max(var)))
        slh.set_v_range_var()
        
        self._viewer.Render()
        
        if filename is not None:
            self._viewer.renwin.save_png(filename)

if __name__ == '__main__':
    dx = 1.
    dy = 1.
    nx = 100
    ny = 100
    Lx = ny * dy
    Ly = nx * dx
    from fipy.meshes.grid2D import Grid2D
    mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)
    from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
    var = DistanceVariable(mesh = mesh, value = -1)
    
    x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
    
    var.setValue(1, where=(x - Lx / 2.)**2 + (y - Ly / 2.)**2 < (Lx / 4.)**2)
    var.calcDistanceFunction()
    viewer = MayaviDistanceViewer(var)
    viewer.plot()
    raw_input("press key to continue")

    var = DistanceVariable(mesh = mesh, value = -1)
    positive = (y > 2. * Ly / 3.) | ((x > Lx / 2.) & (y > Ly / 3.)) | ((y < Ly / 6.) & (x > Lx / 2))
    minY = numerix.minimum.reduce(numerix.compress(positive, y))

    print 'minY',minY
    var.setValue(1, where=positive)
    var.calcDistanceFunction()
    viewer = MayaviDistanceViewer(var)
    viewer.plot()
    raw_input("press key to continue")
