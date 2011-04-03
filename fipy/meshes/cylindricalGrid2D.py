#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "cylindricalGrid2D.py"
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
 #  
 # ###################################################################
 ##

from fipy.tools import parallel

def CylindricalGrid2D(dr=None, dz=None, nr=None, nz=None, dx=1., dy=1., nx=None, ny=None, communicator=parallel):
    from numMesh import cylindricalUniformGrid2D
    from numMesh import cylindricalGrid2D

    from fipy.tools import numerix

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
        return cylindricalUniformGrid2D.CylindricalUniformGrid2D(dx=dx, dy=dy, nx=nx or 1, ny=ny or 1, communicator=communicator)
    else:
        return cylindricalGrid2D.CylindricalGrid2D(dx=dx, dy=dy, nx=nx, ny=ny, communicator=communicator)
