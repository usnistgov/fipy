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

from mayaviViewer import _MayaviViewer

class MayaviVectorViewer(_MayaviViewer):
    
    def __init__(self,vars,title=None,limits={},**kwargs):
        kwargs.update(limits)
        _MayaviViewer.__init__(self,vars,title,**kwargs)
        self.plot()
    
    def _update(self,**kwargs):
        var = kwargs['var']
        src = self.srcs[1]
        src.data.point_data.vectors.to_array()[:]=_MayaviViewer._makeDims3(var.getValue().swapaxes(0,1),move=False)
        src.update()
    
    def _plot(self,**kwargs):
        var = kwargs['var']
        extent = kwargs['extent']
        minDim = kwargs['minDim']
        mesh = var.getMesh()
        dims = mesh.dim
        surf = _MayaviViewer.makeUnstructuredGrid(mesh,minDim/20.)
        from enthought.tvtk.api import tvtk
        surf = _MayaviViewer.makeUnstructuredGrid(mesh,minDim/20.)
        cellCenters = _MayaviViewer._makeDims3(mesh.getCellCenters().swapaxes(0,1),move=False)
        ug = tvtk.UnstructuredGrid(points=cellCenters)
        ug.point_data.vectors = _MayaviViewer._makeDims3(var.getValue().swapaxes(0,1),move=False)
        ug.point_data.vectors.name = 'vectors'
        from enthought.mayavi import mlab
        vecSource = mlab.pipeline.add_dataset(ug)
        gridSource = mlab.pipeline.add_dataset(surf)
        vecs = mlab.pipeline.vectors(vecSource,extent=extent)
        grid = mlab.pipeline.surface(gridSource,extent=extent,opacity=.1)
        self.srcs.extend([gridSource,vecSource])
        self.mods.extend([grid,vecs])

    def _getSuitableVars(self,vars):
        if type(vars) not in [type([]),type(())]:
            vars = [vars]
        vars = [var for var in vars if var.getRank()==1]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError,"Can only plot scalar data"
        return [vars[0]]
