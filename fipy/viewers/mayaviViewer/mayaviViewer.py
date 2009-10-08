#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mayaviViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Daniel Stiles  <daniel.stiles@nist.gov>
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

from mmap import mmap
import os
import subprocess
import tempfile
import time

from fipy.viewers.viewer import _Viewer

class _MayaviViewer(_Viewer):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    The `_MayaviViewer` is the base class for the viewers that use the
    Mayavi_ python plotting package.

    .. Mayavi: http://code.enthought.com/projects/mayavi

    """

    def __init__(self, vars, title=None, **kwlimits):
        """
        Create a `_MayaviViewer`.
        
        :Parameters:
          vars
            a `CellVariable` or tuple of `CellVariable` objects to plot
          title
            displayed at the top of the `Viewer` window
          xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.
        """
        self.vtkdir = tempfile.mkdtemp()
        self.vtkcellfname = os.path.join(self.vtkdir, "cell.vtk")
        self.vtkfacefname = os.path.join(self.vtkdir, "face.vtk")
        self.vtklockfname = os.path.join(self.vtkdir, "lock")

        from fipy.viewers.vtkViewer import VTKCellViewer, VTKFaceViewer

        self.vtkCellViewer = VTKCellViewer(vars=vars)
        self.vtkFaceViewer = VTKFaceViewer(vars=vars)
                                           
#         (self.vtkfile, self.vtkfname) = tempfile.mkstemp('.vtk')
#         (self.mmapfile, self.mmapfname) = tempfile.mkstemp('.mmap')
# #         f2 = os.open(self.mmapfile, os.O_RDWR)
#         os.write(self.mmapfile, '\n' + '\x00' * 1023)
#         self.mmapfile = mmap(self.mmapfile, 1024)
# 
#         from fipy.viewers.vtkViewer import VTKViewer
#         self.vtkViewer = VTKViewer(vars=vars, title=title)

        vars = self.vtkCellViewer.getVars() + self.vtkFaceViewer.getVars()
        
        _Viewer.__init__(self, vars=vars, title=title, **kwlimits)
        
        self.vtkCellViewer.plot(filename=self.vtkcellfname)
        self.vtkFaceViewer.plot(filename=self.vtkfacefname)
        lock = file(self.vtklockfname, 'w')
        lock.close()

        subprocess.Popen(["python", 
                          "/Users/guyer/Documents/research/FiPy/mayavi/fipy/viewers/mayaviViewer/mayaviDaemon.py", 
                          self.vtkcellfname,
                          self.vtkfacefname,
                          self.vtklockfname])

#         self.mods=[]
#         from enthought.mayavi import mlab
#         oldScenes = mlab.get_engine().scenes
#         self.id = len(oldScenes)
# 	if (self.title==''):
#             self.title="FiPy Viewer Window "+str(self.id)
#         self.scene = mlab.figure() #name=self.title)

    def plot(self, filename=None):
        start = time.time()
        plotted = False
        while time.time() - start < 10. and not plotted:
            if not os.path.isfile(self.vtklockfname):
                self.vtkCellViewer.plot(filename=self.vtkcellfname)
                self.vtkFaceViewer.plot(filename=self.vtkfacefname)
                lock = file(self.vtklockfname, 'w')
                if filename is not None:
                    lock.write(filename)
                lock.close()
                plotted = True
        if not plotted:
            print "viewer: NOT READY"
            
#     def plot(self, filename = None):
#         from enthought.mayavi import mlab
# #         mlab.get_engine().current_scene=self.scene
#         xmin = self._getLimit('xmin')
#         xmax = self._getLimit('xmax')
#         ymin = self._getLimit('ymin')
#         ymax = self._getLimit('ymax')
#         zmin = self._getLimit('zmin')
#         zmax = self._getLimit('zmax')
#         datamin = self._getLimit('datamin')
#         datamax = self._getLimit('datamax')
#         for var in self.vars:
#             mesh = var.getMesh()
#             rank = var.getRank()
#             dims = mesh.dim
#             from fipy.tools import numerix
#             x,y,z = None,None,None
#             x = numerix.NUMERIX.min(mesh.getVertexCoords(),axis=1)
#             if dims > 1:
#                 y = x[1:]
#                 x = x[0]
#                 if dims > 2:
#                     z = y[1]
#                     y = y[0]
#             d = None
#             if rank == 0:
#                 d = numerix.NUMERIX.min(var.value)
#             if self._getLimit('xmin') is None and x is not None and (xmin is None or x<xmin):
#                 xmin = x
#             if self._getLimit('ymin') is None and y is not None and (ymin is None or y<ymin):
#                 ymin = y
#             if self._getLimit('zmin') is None and z is not None and (zmin is None or z<zmin):
#                 zmin = z
#             if self._getLimit('datamin') is None and d is not None and (datamin is None or d<datamin):
#                 datamin = d
#             x = numerix.NUMERIX.max(mesh.getVertexCoords(),axis=1)
#             if dims > 1:
#                 y = x[1:]
#                 x = x[0]
#                 if dims > 2:
#                     z = y[1]
#                     y = y[0]
#             d = None
#             if rank == 0:
#                 d = numerix.NUMERIX.max(var.value)
#             if self._getLimit('xmax') is None and x is not None and (xmax is None or x>xmax):
#                 xmax = x
#             if self._getLimit('ymax') is None and y is not None and (ymax is None or y>ymax):
#                 ymax = y
#             if self._getLimit('zmax') is None and z is not None and (zmax is None or z>zmax):
#                 zmax = z
#             if self._getLimit('datamax') is None and d is not None and (datamax is None or d>datamax):
#                 datamax = d
#         if type(xmin) == numerix.ndarray:
#             xmin = xmin[0]
#         if type(xmax) == numerix.ndarray:
#             xmax = xmax[0]
#         if type(ymin) == numerix.ndarray:
#             ymin = ymin[0]
#         if type(ymax) == numerix.ndarray:
#             ymax = ymax[0]
#         if type(zmin) == numerix.ndarray:
#             zmin = zmin[0]
#         if type(zmax) == numerix.ndarray:
#             zmax = zmax[0]
#         if zmax is not None and zmin is not None:
#             minDim = min((xmax-xmin,ymax-ymin,zmax-zmin))
#             maxDim = max((xmax-xmin,ymax-ymin,zmax-zmin))
#         elif ymax is not None and ymin is not None:
#             minDim = min((xmax-xmin,ymax-ymin))
#             maxDim = max((xmax-xmin,ymax-ymin))
#             zmin = 0
#             zmax = minDim/20.
#         else:
#             minDim = xmax-xmin
#             maxDim = xmax-xmin
#             ymin = 0.
#             zmin = 0.
#             ymax = minDim/20.
#             zmax = minDim/20.
#         from fipy.tools.numerix import array
#         for var in self.vars:
#             mesh = var.getMesh()
#             rank = var.getRank()
#             done = False
#             if self.surf is not None:
# 		self._update(var=var,datamin=datamin,datamax=datamax)
#             else:
# 		self._plot(var=var,datamin=datamin,datamax=datamax,extent=[xmin,xmax,ymin,ymax,zmin,zmax],minDim=minDim)
#                 
#         if filename is not None:
#             from enthought.mayavi import mlab
#             mlab.savefig(filename)
    
    def _plot(self):
        pass
    
    def _validFileExtensions(self):
        return [".png",".jpg",".bmp",".tiff",".ps",".eps",".pdf",".rib",".oogl",".iv",".vrml",".obj"]

    def show():
        """Shows this viewer so that the program pauses before displaying the next timestep, or waits to terminate the program till all the windows close"""
        from enthought.mayavi import mlab
        mlab.show()
    show = staticmethod(show)
