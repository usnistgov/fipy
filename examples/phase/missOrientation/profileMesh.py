#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "profileMesh.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:06:52 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##
 
import fipy.tools.dump as dump
from fipy.meshes.grid2D import Grid2D
from fipy.tools.profiler.profiler import calibrate_profiler
from fipy.tools.profiler.profiler import Profiler

def makeMesh(nx, ny):
   return Grid2D(1.,1.,nx,ny)

if __name__ == '__main__':
   fudge = calibrate_profiler(10000)
   profile = Profiler('profile.txt', fudge=fudge)
   mesh = makeMesh(300, 300)
   dump.write(mesh, 'file')
   mesh = dump.read('file')
   profile.stop()
