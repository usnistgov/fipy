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
        self.srcs=[]
    
    _getMul=staticmethod(lambda dims:(dims == 3 and 1 or dims == 2 and 2 or 4))
    
    def _makeDims3(arr,move=True):
        num = arr.shape[0]
        dims = arr.shape[1]
        if dims == 3:
            return arr
        from fipy.tools.numerix import zeros,array
        a = zeros((num*MayaviViewer._getMul(dims),3),t='d')
        a[:,:arr.shape[1]]=array(arr.tolist()*MayaviViewer._getMul(dims))
        if not move:
            return a
        a[num:num*2,2]=1
        if dims == 2:
            return a
        a[num*2:num*3,1]=1
        a[num*3:,2]=1
        a[num*3:,1]=1
        return a
    _makeDims3 = staticmethod(_makeDims3)

    def makeUnstructuredGrid(mesh):
        from fipy.tools.numerix import array
        from enthought.tvtk.api import tvtk
        dims = mesh.dim
        points = mesh.getVertexCoords().swapaxes(0,1)
        numpoints = len(points)
        points = MayaviViewer._makeDims3(points)
        cvi = mesh._getCellVertexIDs().swapaxes(0,1)
        from fipy.tools import numerix
        if (type(cvi)==numerix.ndarray):
            counts = numerix.array([cvi.shape[1]]*cvi.shape[0])
            comp = cvi.flatten()
        else:
            counts = cvi.count(axis=1)
            comp = cvi.compressed()
        lcells = []
        loffsets = []
        total = 0
        mul = MayaviViewer._getMul(dims)
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
        return ug
    makeUnstructuredGrid = staticmethod(makeUnstructuredGrid)
    
    def plot(self, filename = None):
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
            rank = var.getRank()
            dims = mesh.dim
            from fipy.tools import numerix
            x,y,z = None,None,None
            x = numerix.NUMERIX.min(mesh.getVertexCoords(),axis=1)
            if dims > 1:
                y = x[1:]
                x = x[0]
                if dims > 2:
                    z = y[1]
                    y = y[0]
            d = None
            if rank == 0:
                d = numerix.NUMERIX.min(var.value)
            if self._getLimit('xmin') is None and x is not None and (xmin is None or x<xmin):
                xmin = x
            if self._getLimit('ymin') is None and y is not None and (ymin is None or y<ymin):
                ymin = y
            if self._getLimit('zmin') is None and z is not None and (zmin is None or z<zmin):
                zmin = z
            if self._getLimit('datamin') is None and d is not None and (datamin is None or d<datamin):
                datamin = d
            x = numerix.NUMERIX.max(mesh.getVertexCoords(),axis=1)
            if dims > 1:
                y = x[1:]
                x = x[0]
                if dims > 2:
                    z = y[1]
                    y = y[0]
            d = None
            if rank == 0:
                d = numerix.NUMERIX.max(var.value)
            if self._getLimit('xmax') is None and x is not None and (xmax is None or x>xmax):
                xmax = x
            if self._getLimit('ymax') is None and y is not None and (ymax is None or y>ymax):
                ymax = y
            if self._getLimit('zmax') is None and z is not None and (zmax is None or z>zmax):
                zmax = z
            if self._getLimit('datamax') is None and d is not None and (datamax is None or d>datamax):
                datamax = d
        if ymin is None:
            ymin = 0
        if ymax is None:
            ymax = 1
        if zmin is None:
            zmin = 0
        if zmax is None:
            zmax = 1
        from fipy.tools.numerix import array
        for var in self.vars:
            mesh = var.getMesh()
            rank = var.getRank()
            done = False
            if (len(self.srcs)!=0):
                old = self.srcs.pop(0)
                oldMesh = old[0]
                if (oldMesh == mesh):
                    done = True
                    if rank == 0:
                        src = old[1]
                        src.data.cell_data.scalars.to_array()[:] = var.getValue()
                        src.update()
                        s = old[2]
                        s.parent.scalar_lut_manager.data_range=array([datamin,datamax])
                    elif rank == 1:
                        src = old[3]
                        src.data.point_data.vectors.to_array()[:]=MayaviViewer._makeDims3(var.getValue().swapaxes(0,1),move=False)
                        src.update()
                    self.srcs.append(old)
            if (not done):
                dims = mesh.dim
                from enthought.tvtk.api import tvtk
                surf = MayaviViewer.makeUnstructuredGrid(mesh)
                if rank == 0:
                    surf.cell_data.scalars = var.getValue()
                    surf.cell_data.scalars.name = 'scalars'                    
                    from enthought.mayavi import mlab
                    src = mlab.pipeline.add_dataset(surf)
                    s = mlab.pipeline.surface(src,extent=[xmin, xmax, ymin, ymax, zmin, zmax],vmin=datamin,vmax=datamax)
                    self.srcs.append([mesh,src,s])
                elif rank == 1:
                    cellCenters = MayaviViewer._makeDims3(mesh.getCellCenters().swapaxes(0,1),move=False)
                    ug = tvtk.UnstructuredGrid(points=cellCenters)
                    ug.point_data.vectors = MayaviViewer._makeDims3(var.getValue().swapaxes(0,1),move=False)
                    ug.point_data.vectors.name = 'vectors'
                    from enthought.mayavi import mlab
                    vecSource = mlab.pipeline.add_dataset(ug)
                    gridSource = mlab.pipeline.add_dataset(surf)
                    vecs = mlab.pipeline.vectors(vecSource,extent=[xmin, xmax, ymin, ymax, zmin, zmax] )
                    grid = mlab.pipeline.surface(gridSource,extent=[xmin,xmax,ymin,ymax,zmin,zmax],opacity=.1)
                    self.srcs.append([mesh,gridSource,grid,vecSource,vecs])
                else:
                    raise TypeError("The data for mayavi must be scalars or vectors")
                
        if filename is not None:
            from enthought.mayavi import mlab
            mlab.savefig(filename)
        def _validFileExtensions(self):
            return [".png",".jpg",".bmp",".tiff",".ps",".eps",".pdf",".rib",".oogl",".iv",".vrml",".obj"]

    def show(self):
        from enthought.mayavi import mlab
        mlab.show()
