#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mayaviSurfactantViewer.py"
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

from fipy.viewers.viewer import AbstractViewer
from fipy.viewers import MeshDimensionError
from fipy.tools import numerix

class MayaviSurfactantViewer(AbstractViewer):

    """
    The `MayaviSurfactantViewer` creates a viewer with the Mayavi_ python
    plotting package that displays a `DistanceVariable`.

    .. _Mayavi: http://mayavi.sourceforge.net/

    """

    def __init__(self, distanceVar, surfactantVar=None, levelSetValue=0., title=None, smooth=0, zoomFactor=1., animate=False, limits={}, **kwlimits):
        """
        Create a `MayaviSurfactantViewer`.

            >>> from fipy import *
            >>> dx = 1.
            >>> dy = 1.
            >>> nx = 11
            >>> ny = 11
            >>> Lx = ny * dy
            >>> Ly = nx * dx
            >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)
            >>> # from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
            >>> var = DistanceVariable(mesh = mesh, value = -1)

            >>> x, y = mesh.cellCenters

            >>> var.setValue(1, where=(x - Lx / 2.)**2 + (y - Ly / 2.)**2 < (Lx / 4.)**2)
            >>> var.calcDistanceFunction()
            >>> viewer = MayaviSurfactantViewer(var, smooth = 2)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

            >>> var = DistanceVariable(mesh = mesh, value = -1)

            >>> var.setValue(1, where=(y > 2. * Ly / 3.) | ((x > Lx / 2.) & (y > Ly / 3.)) | ((y < Ly / 6.) & (x > Lx / 2)))
            >>> var.calcDistanceFunction()
            >>> viewer = MayaviSurfactantViewer(var)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

            >>> viewer = MayaviSurfactantViewer(var, smooth = 2)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

        :Parameters:

          - `distanceVar`: a `DistanceVariable` object.
          - `levelSetValue`: the value of the contour to be displayed
          - `title`: displayed at the top of the `Viewer` window
          - `animate`: whether to show only the initial condition and the
          - `limits`: a dictionary with possible keys `xmin`, `xmax`,
            `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.  A 1D
            `Viewer` will only use `xmin` and `xmax`, a 2D viewer will also
            use `ymin` and `ymax`, and so on.  All viewers will use
            `datamin` and `datamax`.  Any limit set to a (default) value of
            `None` will autoscale.
            moving top boundary or to show all contours (Default)
        """

        kwlimits.update(limits)
        AbstractViewer.__init__(self, vars=[], title=title, **kwlimits)
        import mayavi
        self._viewer = mayavi.mayavi()
        self.distanceVar = distanceVar
        if surfactantVar is None:
            self.surfactantVar = numerix.zeros(len(self.distanceVar), 'd')
        else:
            self.surfactantVar = surfactantVar
        self.smooth = smooth
        self.zoomFactor = zoomFactor

        self.animate = animate
        if animate:
            self._initialCondition = None

        if distanceVar.mesh.dim != 2:
            raise MeshDimensionError('The MayaviIsoViewer only works for 2D meshes.')

    def _getStructure(self):

        ##maxX = self.distanceVar.mesh.faceCenters[0].max()
        ##minX = self.distanceVar.mesh.faceCenters[0].min()

        IDs = numerix.nonzero(self.distanceVar._cellInterfaceFlag)[0]
        coordinates = numerix.take(numerix.array(self.distanceVar.mesh.cellCenters).swapaxes(0,1), IDs)

        coordinates -= numerix.take(numerix.array(self.distanceVar.grad * self.distanceVar).swapaxes(0,1), IDs)

        coordinates *= self.zoomFactor

        shiftedCoords = coordinates.copy()
        shiftedCoords[:,0] = -coordinates[:,0] ##+ (maxX - minX)
        coordinates = numerix.concatenate((coordinates, shiftedCoords))

        from lines import _getOrderedLines

        lines = _getOrderedLines(range(2 * len(IDs)), coordinates, thresholdDistance = self.distanceVar.mesh._cellDistances.min() * 10)

        data = numerix.take(self.surfactantVar, IDs)

        data = numerix.concatenate((data, data))

        tmpIDs = numerix.nonzero(data > 0.0001)[0]
        if len(tmpIDs) > 0:
            val = numerix.take(data, tmpIDs).min()
        else:
            val = 0.0001

        data = numerix.where(data < 0.0001,
                             val,
                             data)

        for line in lines:
            if len(line) > 2:
                for smooth in range(self.smooth):
                    for arr in (coordinates, data):
                        tmp = numerix.take(arr, line)
                        tmp[1:-1] = tmp[2:] * 0.25 + tmp[:-2] * 0.25 + tmp[1:-1] * 0.5
                        if len(arr.shape) > 1:
                            for i in range(len(arr[0])):
                                arrI = arr[:,i].copy()
                                numerix.put(arrI, line, tmp[:,i])
                                arr[:,i] = arrI
                        else:
                            numerix.put(arrI, line, tmp)

        name = self.title
        name = name.strip()
        if name == '':
            name = None

        coords = numerix.zeros((coordinates.shape[0], 3), 'd')
        coords[:,:coordinates.shape[1]] = coordinates

        import pyvtk

        ## making lists as pyvtk doesn't know what to do with numpy arrays

        coords = list(coords)
        coords = map(lambda coord: [float(coord[0]),float(coord[1]), float(coord[2])], coords)

        data = list(data)
        data = map(lambda item: float(item), data)

        return (pyvtk.UnstructuredGrid(points = coords,
                                       poly_line = lines),
                pyvtk.PointData(pyvtk.Scalars(data, name = name)))

    def plot(self, filename = None):

        structure, data = self._getStructure()

        import pyvtk
        import tempfile
        import os

        if self.animate:
            if self._initialCondition is None:
                data = pyvtk.VtkData(structure, 0)
                (f, tempFileName) = tempfile.mkstemp('.vtk')
                data.tofile(tempFileName)
                self._viewer.open_vtk(tempFileName, config=0)
                os.close(f)
                os.remove(tempFileName)

                self._viewer.load_module('SurfaceMap', 0)

                rw = self._viewer.get_render_window()
                rw.z_plus_view()

                self._initialCondition = self._viewer.get_current_dvm_name()
            else:
                self._viewer.mayavi.del_dvm(self._viewer.get_current_dvm_name())

        data = pyvtk.VtkData(structure, data)

        (f, tempFileName) = tempfile.mkstemp('.vtk')
        data.tofile(tempFileName)
        self._viewer.open_vtk(tempFileName, config=0)

        os.close(f)
        os.remove(tempFileName)
        self._viewer.load_module('SurfaceMap', 0)

        if not self.animate:
            rw = self._viewer.get_render_window()
            rw.z_plus_view()

        ## display legend
        dvm = self._viewer.get_current_dvm()
        mm = dvm.get_current_module_mgr()
        slh = mm.get_scalar_lut_handler()
        slh.legend_on.set(1)
        slh.legend_on_off()

        ## display legend with correct range
        slh.range_on_var.set(1)
        slh.v_range_on_var.set(1)

        xmax = self._getLimit('datamax', default=self.surfactantVar.max())
        xmin = self._getLimit('datamin', default=self.surfactantVar.min())

        slh.range_var.set((xmin, xmax))
        slh.set_range_var()

        slh.v_range_var.set((float(self.surfactantVar.min()), float(self.surfactantVar.max())))
        slh.set_v_range_var()

        self._viewer.Render()

        if filename is not None:
            self._viewer.renwin.save_png(filename)

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
