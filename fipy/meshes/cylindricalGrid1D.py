#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "cylindricalGrid1D.py"
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

def CylindricalGrid1D(dr=None, nr=None, dx=1., nx=None):
    from numMesh import cylindricalUniformGrid1D
    from numMesh import cylindricalGrid1D

    from fipy.tools import numerix

    if dr is not None:
        dx = dr
        
    nx = nr or nx

    if numerix.getShape(dx) == ():
        return cylindricalUniformGrid1D.CylindricalUniformGrid1D(dx=dx, nx=nx or 1)
    else:
        return cylindricalGrid1D.CylindricalGrid1D(dx=dx, nx=nx)



