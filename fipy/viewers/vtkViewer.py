#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vtkCellViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Stiles  <daniel.stiles@nist.gov>
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

import os
import time

from fipy.viewers.viewer import _Viewer

class VTKCellViewer(_Viewer):
    """Renders `CellVariable` data in VTK format
    """
    def __init__(self, vars, title=None, limits={}, **kwlimits):
        """Creates a VTKCellViewer

        :Parameters:
          vars
            a `CellVariable` or a tuple of them
          title
            displayed at the top of the `Viewer` window
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to 
            a (default) value of `None` will autoscale.
        """
        kwlimits.update(limits)
        _Viewer.__init__(self, vars=vars, title=title, **kwlimits)

        mesh = self.vars[0].getMesh()
        
        (self.vtkDataSet,
         self.vtkFaceOffset,
         self.vtkCellOffset,
         self.vtkTotalPoints) = mesh.getVTKDataSet()

        for var in self.vars:
            name = self._getArrayName(var)

            rank = var.getRank()
            if rank == 0:
                value = var.getValue()
                i = self.vtkDataSet.cell_data.add_array(value)
                self.vtkDataSet.cell_data.get_array(i).name = name
                self.vtkDataSet.cell_data.set_active_scalars(name)
            else:
                value = numerix.zeros((self.vtkTotalPoints,) + rank * (3,), 'd')
                value[self.vtkCellOffset:] = var.getValue().swapaxes(-2,-1)
            
                i = self.vtkDataSet.point_data.add_array(value)
                self.vtkDataSet.point_data.get_array(i).name = name

                if rank == 0:
                    self.vtkDataSet.point_data.set_active_scalars(name)
                elif rank == 1:
                    self.vtkDataSet.point_data.set_active_vectors(name)
                else:
                    self.vtkDataSet.point_data.set_active_tensors(name)
        
    @staticmethod
    def _getArrayName(var):
        return var.name or "%s #%d" % (var.__class__.__name__, id(var))
    
    def plot(self, filename=None):
        for var in self.vars:
            name = self._getArrayName(var)
            rank = var.getRank()
            if rank == 0:
                value = var.getValue()
                self.vtkDataSet.cell_data.get_array(name).to_array()[:] = value
            else:
                value = numerix.zeros((self.vtkTotalPoints,) + rank * (3,), 'd')
                value[self.vtkCellOffset:] = var.getValue().swapaxes(-2,-1)

                self.vtkDataSet.point_data.get_array(name).to_array()[:] = value

        from enthought.tvtk.api import tvtk
        w = tvtk.UnstructuredGridWriter(input=self.vtkDataSet, file_name=filename)
        w.write()

#         w = tvtk.XMLUnstructuredGridWriter(input=self.vtkDataSet, file_name=filename) #[:-1] + "u")
#         w.write()

#         from enthought.tvtk.misc import write_data
#         write_data(self.vtkDataSet, filename)

    def _getSuitableVars(self,vars):
        if type(vars) not in [type([]),type(())]:
            vars = [vars]
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in vars if isinstance(var, CellVariable)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError,"Can only plot CellVariable data"
        vars = [var for var in vars if var.getMesh()==vars[0].getMesh()]
        return vars

if __name__ == "__main__": 
#     import fipy.tests.doctestPlus
#     fipy.tests.doctestPlus.execButNoTest()

    from fipy import *
    m = Grid3D(nx=3, ny=4, nz=5)
    x, y, z = m.getCellCenters()
    v1 = CellVariable(mesh=m, value=x*y*z, name="v1")
    v2 = CellVariable(mesh=m, value=x*y*y, name="v2")
    v3 = v1.getGrad()
    v3.name = "v3"
#     vw = VTKCellViewer(vars=(v1, v2))
    vw = VTKCellViewer(vars=(v1, v2, v3))
    
    vw.plot(filename="vtk.vtk")

#     m = Grid2D(nx=1, ny=2)
#     x, y = m.getCellCenters()
#     v1 = CellVariable(mesh=m, value=x*y, name="v1")
#     v2 = CellVariable(mesh=m, value=x*x) #, name="v2")
#     vw = VTKCellViewer(vars=(v1, v2))

#     m = Grid1D(nx=10)
#     x,  = m.getCellCenters()
#     v1 = CellVariable(mesh=m, value=x*x, name="v1")
#     v2 = CellVariable(mesh=m, value=x) #, name="v2")
#     vw = VTKCellViewer(vars=(v1, v2))

