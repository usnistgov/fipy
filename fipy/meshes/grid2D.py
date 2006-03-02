#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "grid2D.py"
 #                                    created: 11/20/03 {4:47:54 PM} 
 #                                last update: 3/2/06 {3:46:32 PM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-20 JEG 1.0 original
 # ###################################################################
 ##

from numMesh import uniformGrid2D
from numMesh import grid2D

from fipy.tools import numerix

def Grid2D(dx = 1., dy = 1., nx = None, ny = None):
    if numerix.getShape(dx) == () and numerix.getShape(dy) == ():
        return uniformGrid2D.UniformGrid2D(dx = dx, dy = dy, 
                                           nx = nx or 1, ny = ny or 1)
    else:
        return grid2D.Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)
