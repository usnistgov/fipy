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

class _MayaviViewer(_Viewer):
        
    def __init__(self, vars, title=None, **kwlimits):
        _Viewer.__init__(self, vars=vars, title=title, **kwlimits)
        self.srcs=[]
        self.mods=[]
        from enthought.mayavi import mlab
        oldScenes = mlab.get_engine().scenes
        self.id = len(oldScenes)
	if (self.title==''):
            self.title="FiPy Viewer Window "+str(self.id)
        self.scene = mlab.figure(name=self.title)

    def createColorbar(self):
        from enthought.mayavi import mlab
        mlab.colorbar(object=self.mods[0],orientation='horizontal',nb_labels=5)
        
    _getMul=staticmethod(lambda dims,mul=True:((not mul or dims == 3) and 1 or dims == 2 and 2 or 4))
    
    def _makeDims3(arr,expand=1.,move=True):
        '''This method takes a numpy array and expands it from some dimension to 3 dimensions.
        :Parameters:
          arr
            The array to expand
          expand
            How far to move the copies of the original points (only used if move is True
          move
            Whether or not to copy the points and translate them.  This will make 4 copies of a 1D array and 2 copies of a 2D array.'''
        num = arr.shape[0]
        dims = arr.shape[1]
        if dims == 3:
            return arr
        from fipy.tools.numerix import zeros,array
        a = zeros((num*_MayaviViewer._getMul(dims,move),3),t='d')
        a[:,:arr.shape[1]]=array(arr.tolist()*_MayaviViewer._getMul(dims,move))
        if not move:
            return a
        a[num:num*2,2]=expand
        if dims == 2:
            return a
        a[num*2:num*3,1]=expand
        a[num*3:,2]=expand
        a[num*3:,1]=expand
        return a
    _makeDims3 = staticmethod(_makeDims3)

    def makeUnstructuredGrid(mesh,minDim):
        '''This method takes a fipy mesh and turns it into an enthought.tvtk unstructured grid.
        :Parameters:
          mesh
            the mesh to make into an unstructured gird
          minDim
            tells the method how far to extrude out 1D and 2D meshes'''
        from fipy.tools.numerix import array
        from enthought.tvtk.api import tvtk
        dims = mesh.dim
        points = mesh.getVertexCoords().swapaxes(0,1)
        numpoints = len(points)
        points = _MayaviViewer._makeDims3(points,expand=minDim)
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
        mul = _MayaviViewer._getMul(dims)
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
        from enthought.mayavi import mlab
        mlab.get_engine().current_scene=self.scene
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
        if type(xmin) == numerix.ndarray:
            xmin = xmin[0]
        if type(xmax) == numerix.ndarray:
            xmax = xmax[0]
        if type(ymin) == numerix.ndarray:
            ymin = ymin[0]
        if type(ymax) == numerix.ndarray:
            ymax = ymax[0]
        if type(zmin) == numerix.ndarray:
            zmin = zmin[0]
        if type(zmax) == numerix.ndarray:
            zmax = zmax[0]
        if zmax is not None and zmin is not None:
            minDim = min((xmax-xmin,ymax-ymin,zmax-zmin))
            maxDim = max((xmax-xmin,ymax-ymin,zmax-zmin))
        elif ymax is not None and ymin is not None:
            minDim = min((xmax-xmin,ymax-ymin))
            maxDim = max((xmax-xmin,ymax-ymin))
            zmin = 0
            zmax = minDim/20.
        else:
            minDim = xmax-xmin
            maxDim = xmax-xmin
            ymin = 0.
            zmin = 0.
            ymax = minDim/20.
            zmax = minDim/20.
        from fipy.tools.numerix import array
        for var in self.vars:
            mesh = var.getMesh()
            rank = var.getRank()
            done = False
            if (len(self.srcs)!=0):
		self._update(var=var,datamin=datamin,datamax=datamax)
            else:
		self._plot(var=var,datamin=datamin,datamax=datamax,extent=[xmin,xmax,ymin,ymax,zmin,zmax],minDim=minDim)
                
        if filename is not None:
            from enthought.mayavi import mlab
            mlab.savefig(filename)
    
    def _plot(self):
        pass
    
    def _validFileExtensions(self):
        return [".png",".jpg",".bmp",".tiff",".ps",".eps",".pdf",".rib",".oogl",".iv",".vrml",".obj"]

    def show():
        from enthought.mayavi import mlab
        mlab.show()
    show = staticmethod(show)
