#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gistViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 10/25/06 {4:15:28 PM} { 2:45:36 PM}
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

from fipy.viewers.gistViewer.gistViewer import GistViewer

from fipy.variables.cellVariable import CellVariable
from fipy.variables.faceVariable import FaceVariable

from fipy.tools import numerix

class GistVectorViewer(GistViewer):
    
    def __init__(self, vars, title = ''):
        """
            >>> from fipy import *
            >>> mesh = Grid2D(nx=50, ny=100, dx=0.1, dy=0.01)
            >>> x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
            >>> var = CellVariable(mesh=mesh, name=r"$sin(x y)$", value=numerix.sin(x * y))
            >>> viewer = GistVectorViewer(vars=var.getGrad(), 
            ...                           limits={'ymin':0.1, 'ymax':0.9, 'datamin':-0.9, 'datamax':2.0},
            ...                           title="GistVectorViewer test")
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

            >>> viewer = GistVectorViewer(vars=var.getFaceGrad(), 
            ...                           limits={'ymin':0.1, 'ymax':0.9, 'datamin':-0.9, 'datamax':2.0},
            ...                           title="GistVectorViewer test")
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer
            
            >>> mesh = Tri2D(nx=50, ny=100, dx=0.1, dy=0.01)
            >>> x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
            >>> var = CellVariable(mesh=mesh, name=r"$sin(x y)$", value=numerix.sin(x * y))
            >>> viewer = GistVectorViewer(vars=var.getGrad(), 
            ...                           limits={'ymin':0.1, 'ymax':0.9, 'datamin':-0.9, 'datamax':2.0},
            ...                           title="GistVectorViewer test")
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

            >>> viewer = GistVectorViewer(vars=var.getFaceGrad(), 
            ...                           limits={'ymin':0.1, 'ymax':0.9, 'datamin':-0.9, 'datamax':2.0},
            ...                           title="GistVectorViewer test")
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

        """
	GistViewer.__init__(self, vars=vars, title=title)
        
    def _getSuitableVars(self, vars):
        vars = [var for var in GistViewer._getSuitableVars(self, vars) \
          if (var.getMesh().getDim() == 2 \
              and (isinstance(var, FaceVariable) \
                   or isinstance(var, CellVariable)) and var.getRank() == 1)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "Can only plot 2D vector data"
        # this viewer can only display one variable
        return [vars[0]]
        
    def plot(self, filename = None):
        import gist

        gist.window(self.id, wait = 1)
        gist.pltitle(self.title)
        gist.animate(1)
        
        var = self.vars[0]
        
        if isinstance(var, FaceVariable):
            x, y = var.getMesh().getFaceCenters()
        elif isinstance(var, CellVariable):
            x, y = var.getMesh().getCellCenters()
        
        gist.plmesh(numerix.array([y, y]), numerix.array([x, y]))

        vx = numerix.array(var[0])
        vy = numerix.array(var[1])
        
        maxVec = var.getMag().max()
        maxGrid = var.getMesh()._getCellDistances().max()
        
        gist.plv(numerix.array([vy,vy]), numerix.array([vx,vx]), scale=maxGrid / maxVec * 3, hollow=1, aspect=0.25) #,scale=0.002)
        
        if filename is not None:
            
            gist.hcp_file(filename)
            gist.hcp()

        gist.fma()
        
    def getArray(self):
        pass
        
if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
