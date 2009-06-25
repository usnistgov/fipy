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
        from enthought.mayavi.core.module_manager import ModuleManager
        self.modMan = ModuleManager()
        self.modMan.scalar_lut_manager.use_default_range=False
        
    _getMul=lambda self,dims:(dims == 3 and 1 or dims == 2 and 2 or 4)
    
    def _makeDims3(self, arr,move=True):
        num = arr.shape[0]
        dims = arr.shape[1]
        if dims == 3:
            return arr
        from fipy.tools.numerix import zeros
        a = zeros((num*self._getMul(dims),3))
        a[:,:arr.shape[1]]=list(arr)*self._getMul(dims)
        if not move:
            return a
        a[num:num*2,2]=1
        if dims == 2:
            return a
        a[num*2:num*3,1]=1
        a[num*3:,2]=1
        a[num*3:,1]=1
        return a
    
    def plot(self, filename = None):
        if not self.e.scenes:
            self.e.new_scene()
        xmin = self._getLimit('xmin')
        xmax = self._getLimit('xmax')
        ymin = self._getLimit('ymin')
        ymax = self._getLimit('ymax')
        zmin = self._getLimit('zmin')
        zmax = self._getLimit('zmax')
        datamin = self._getLimit('datamin')
        datamax = self._getLimit('datamax')
        for var in self.vars:
            mesh = var.getMesh()
            from fipy.tools import numerix
            x,y,z = numerix.NUMERIX.min(mesh.getVertexCoords(),axis=1)
            d = numerix.NUMERIX.min(var.value)
            if self._getLimit('xmin') is None and (xmin is None or x<xmin):
                xmin = x
            if self._getLimit('ymin') is None and (ymin is None or y<ymin):
                ymax = y
            if self._getLimit('zmin') is None and (zmin is None or z<zmin):
                zmax = z
            if self._getLimit('datamin') is None and (datamin is None or d<datamin):
                datamin = d
            x,y,z = numerix.NUMERIX.max(mesh.getVertexCoords(),axis=1)
            d = numerix.NUMERIX.max(var.value)
            if self._getLimit('xmax') is None and (xmax is None or x>xmax):
                xmax = x
            if self._getLimit('ymax') is None and (ymax is None or y>ymax):
                ymax = y
            if self._getLimit('zmax') is None and (zmax is None or z>zmax):
                zmax = z
            if self._getLimit('datamax') is None and (datamax is None or d>datamax):
                datamax = d
        self.modMan.scalar_lut_manager.data_range = (datamin,datamax)
        from fipy.tools.numerix import array
        for var in self.vars:
            mesh = var.getMesh()
            name = var.getName()
            if name is '':
                name = 'default'

            rank = var.getRank()
            dims = mesh.dim
            ug = None
            mod = None
            from enthought.tvtk.api import tvtk
            if rank == 0:
                points = mesh.getVertexCoords().swapaxes(0,1)
                numpoints = len(points)
                points = self._makeDims3(points)
                cvi = mesh._getCellVertexIDs().swapaxes(0,1)
                from fipy.tools import numerix
                if (type(cvi)==numerix.ma.masked_array):
                    counts = cvi.count(axis=1)
                    comp = cvi.compressed()
                else:
                    counts = [cvi.shape[1]]*cvi.shape[0]
                    comp = cvi.reshape(-1)
                lcells = []
                loffsets = []
                total = 0
                mul = self._getMul(dims)
                num = 0
                for n in counts:
                    loffsets += [total*mul+num]
                    lcells += [n*mul]
                    for i in range(n):
                        lcells += [comp[total+i]]
                        if (dims < 3):
                            lcells += [comp[total+i]+numpoints]
                            if (dims < 2):
                                lcells += [comp[total+i]+numpoints*2]
                                lcells += [comp[total+i]+numpoints*3]
                    total += n
                    num += 1
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
                mod.module_manager = self.modMan
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
            from enthought.mayavi import mlab
            mlab.savefig(filename)
        def _validFileExtensions(self):
            return [".png",".jpg",".bmp",".tiff",".ps",".eps",".pdf",".rib",".oogl",".iv",".vrml",".obj"]

    def show(self):
        from enthought.mayavi import mlab
        mlab.show()
