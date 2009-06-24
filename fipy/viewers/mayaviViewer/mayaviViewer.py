#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "matplotlibViewer.py"
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

from fipy.viewers.viewer import _Viewer

class MayaviViewer(_Viewer):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    The `_MayaviViewer` is the base class for the viewers that use the
    Mayavi python plotting package.


    """
        
    def __init__(self, vars, title=None, **kwlimits):
        """
        Create a `MayaviViewer`.
        
        :Parameters:
          vars
            a `CellVariable` or tuple of `CellVariable` objects to plot
          title
            displayed at the top of the `Viewer` window
          figaspect
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is an array, figaspect will determine the width and
            height for a figure that would fit array preserving aspect ratio.
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.
        """
        _Viewer.__init__(self, vars=vars, title=title, **kwlimits)
        self.oldSrcs = []
        self.oldMods = []
        from enthought.mayavi.api import Engine
        self.e = Engine()
        self.e.start()
        self.e.new_scene(title)

    def plot(self, filename = None):
        from fipy.tools.numerix import array
        for var in self.vars:
            mesh = var.getMesh()
            name = var.getName()
            if name is '':
                name = 'default'

            rank = var.getRank()
            ug = None
            mod = None
            if rank == 0:
                from enthought.tvtk.api import tvtk
                cvi = mesh._getCellVertexIDs().swapaxes(0,1)
                counts = cvi.count(axis=1)
                comp = cvi.compressed()
                lcells = []
                loffsets = []
                total = 0
                num = 0
                for n in counts:
                    loffsets += [total+num]
                    lcells += [n]
                    for i in range(n):
                        lcells += [comp[total+i]]
                    total += n
                    num += 1
                points = mesh.getVertexCoords().swapaxes(0,1)
                cells = array(lcells)
                offset = array(loffsets)
                cps_type = tvtk.ConvexPointSet().cell_type
                cell_types = array([cps_type]*num)
                cell_array = tvtk.CellArray()
                cell_array.set_cells(num, cells)
                
                ug = tvtk.UnstructuredGrid(points=points)
                ug.set_cells(cell_types, offset, cell_array)
                ug.cell_data.scalars = var.value
                ug.cell_data.scalars.name = 'scalars'

                from enthought.mayavi.modules.api import Surface
                mod = Surface()
                
            elif rank == 1:
                cellCenters = mesh.getCellCenters().swapaxes(0,1)
                ug = tvtk.UnstructuredGrid(points=cellCenters)
                ug.point_data.vectors = var.value
                ug.point_data.vectors.name = 'vectors'
                
                from enthought.mayavi.modules import Vectors
                mod = Vectors()

            if (mod is None) or (ug is None):
                raise TypeError("The data for mayavi must be scalars or vectors")

            from enthought.mayavi.sources.api import VTKDataSource
            src = VTKDataSource(data=ug)
            self.e.add_source(src)
            self.e.add_module(mod,obj=src)
            if (len(self.oldSrcs) > 0):
                oldSrc = self.oldSrcs.pop(0)
                oldMod = self.oldMods.pop(0)
                oldMod.stop()
                oldSrc.stop()
            self.oldSrcs += [src]
            self.oldMods += [mod]

        if filename is not None:
            pass
        def _validFileExtensions(self):
            return [".png"]

    def show(self):
        from enthought.mayavi import mlab
        mlab.show()
