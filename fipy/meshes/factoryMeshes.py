#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Alexander Mont <alexander.mont@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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
 
from fipy.tools import parallel
from fipy.tools import numerix

__all__ = ["Grid3D", "Grid2D", "Grid1D", "CylindricalGrid2D", "CylindricalGrid1D"]

def Grid3D(dx = 1., dy = 1., dz = 1., nx = None, ny = None, nz = None, overlap=2, communicator=parallel):
    import uniformGrid3D
    import grid3D

    if numerix.getShape(dx) == () \
      and numerix.getShape(dy) == () \
      and numerix.getShape(dz) == ():
        if nx is None:
            nx = 1
        if ny is None:
            ny = 1
        if nz is None:
            nz = 1
        return uniformGrid3D.UniformGrid3D(dx = dx, dy = dy, dz = dz,
                                           nx = nx or 1, ny = ny or 1, nz = nz or 1,
                                           overlap=overlap, communicator=communicator)
    else:
        return grid3D.Grid3D(dx = dx, dy = dy, dz = dz, nx = nx, ny = ny, nz = nz,
                             overlap=overlap, communicator=communicator) 

def Grid2D(dx=1., dy=1., nx=None, ny=None, overlap=2, communicator=parallel):
    import uniformGrid2D
    import grid2D

    if numerix.getShape(dx) == () and numerix.getShape(dy) == ():
        if nx is None:
            nx = 1
        if ny is None:
            ny = 1
        return uniformGrid2D.UniformGrid2D(dx=dx, dy=dy, 
                                           nx=nx, ny=ny,
                                           overlap=overlap,
                                           communicator=communicator)
    else:
        return grid2D.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny, overlap=overlap, communicator=communicator)

def Grid1D(dx=1., nx=None, overlap=2, communicator=parallel):
    import uniformGrid1D
    import grid1D
    
    if numerix.getShape(dx) == ():
        if nx is None:
            nx = 1
        return uniformGrid1D.UniformGrid1D(dx=dx, nx=nx, overlap=overlap, communicator=communicator)
    else:
        return grid1D.Grid1D(dx=dx, nx=nx, overlap=overlap, communicator=communicator)

def CylindricalGrid2D(dr=None, dz=None, 
                      nr=None, nz=None, 
                      dx=1., dy=1., 
                      nx=None, ny=None,
                      origin=((0,),(0,)),
                      overlap=2,
                      communicator=parallel):
    import cylindricalUniformGrid2D
    import cylindricalGrid2D

    if dr is not None:
        dx = dr

    if dz is not None:
        dy = dz

    nx = nr or nx
    ny = nz or ny

    if dx is None:
        dx = 1.
    if dy is None:
        dy = 1.
    
    if numerix.getShape(dx) == () and numerix.getShape(dy) == ():
        return cylindricalUniformGrid2D.CylindricalUniformGrid2D(dx=dx, dy=dy, nx=nx or 1, ny=ny or 1, origin=origin, overlap=overlap, communicator=communicator)
    else:
        return cylindricalGrid2D.CylindricalGrid2D(dx=dx, dy=dy, nx=nx, ny=ny, origin=origin, overlap=overlap, communicator=communicator)

def CylindricalGrid1D(dr=None, nr=None, dx=1., nx=None, origin=(0,), overlap=2, communicator=parallel):
    import cylindricalUniformGrid1D
    import cylindricalGrid1D

    if dr is not None:
        dx = dr
        
    nx = nr or nx

    if numerix.getShape(dx) == ():
        return cylindricalUniformGrid1D.CylindricalUniformGrid1D(dx=dx, nx=nx or 1, origin=origin, overlap=overlap, communicator=parallel)
    else:
        return cylindricalGrid1D.CylindricalGrid1D(dx=dx, nx=nx, origin=origin, overlap=overlap, communicator=parallel)

