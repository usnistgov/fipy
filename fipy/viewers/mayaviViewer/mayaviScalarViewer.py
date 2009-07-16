#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mayaviScalarViewer.py"
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

from mayaviViewer import _MayaviViewer

class MayaviScalarViewer(_MayaviViewer):
    
    def __init__(self,vars,title=None,limits={},**kwlimits):
        kwlimits.update(limits)
        _MayaviViewer.__init__(self,vars,title,**kwlimits)
        self.plot()

    def _update(self,**kwargs):
        var = kwargs['var']
        datamin=kwargs['datamin']
        datamax=kwargs['datamax']
        src = self.srcs[0]
        src.data.cell_data.scalars.to_array()[:]=var.getValue()
        src.update()
        s = self.mods[0]
	from fipy.tools.numerix import array
        s.parent.scalar_lut_manager.data_range=array([datamin,datamax])
    
    def _plot(self,**kwargs):
        var = kwargs['var']
        extent = kwargs['extent']
        minDim = kwargs['minDim']
        datamin = kwargs['datamin']
        datamax = kwargs['datamax']
        mesh = var.getMesh()
        dims = mesh.dim
        surf = _MayaviViewer.makeUnstructuredGrid(mesh,minDim/20.)
        surf.cell_data.scalars = var.getValue()
        surf.cell_data.scalars.name = 'scalars'                    
        from enthought.mayavi import mlab
        src = mlab.pipeline.add_dataset(surf)
        s = mlab.pipeline.surface(src,extent=extent,vmin=datamin,vmax=datamax)
        self.srcs.append(src)
        self.mods.append(s)
        self.createColorbar()

    def _getSuitableVars(self,vars):
        if type(vars) not in [type([]),type(())]:
            vars = [vars]
        vars = [var for var in vars if var.getRank()==0]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError,"Can only plot scalar data"
        return [vars[0]]
