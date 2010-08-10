#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "grid2D.py"
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

def Grid2D(dx=1., dy=1., nx=None, ny=None, overlap=2, communicator=parallel):
    from numMesh import uniformGrid2D
    from numMesh import grid2D

    from fipy.tools import numerix
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
