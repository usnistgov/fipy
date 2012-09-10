#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "levelset.py"
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

__docformat__ = 'restructuredtext'

r"""

Wrapper functions to provide a common interface for lsmlib and
scikit-fmm.

"""

from fipy.tests.doctestPlus import register_skipper
from fipy.tools import numerix
import sys
import os

def _parseLSMSolver():
    args = [s.lower() for s in sys.argv[1:]]
    # any command-line specified solver takes precedence over environment variables
    if '--lsmlib' in args:
        return "lsmlib"
    elif '--skfmm' in args:
        return "skfmm"
    elif 'FIPY_LSM' in os.environ:
        return os.environ['FIPY_LSM'].lower()
    else:
        return None

LSM_SOLVER = _parseLSMSolver()

if LSM_SOLVER is None:
    try:
        import pylsmlib
        LSM_SOLVER = 'lsmlib'
    except Exception:    
        try:
            import skfmm
            LSM_SOLVER = 'skfmm'
        except Exception:
            pass

def _checkForLSMLIB():
    return LSM_SOLVER == 'lsmlib'

def _checkForLSM():
    return LSM_SOLVER != None

register_skipper(flag="LSM",
                 test=_checkForLSM,
                 why="Neither `lsmlib` nor `skfmm` can be found on the $PATH")

register_skipper(flag="LSMORDER1",
                 test=_checkForLSMLIB,
                 why="only `lsmlib` can perform first order level set calculations")
                     
def calcDistanceFunction(phi, mesh, order):
    
    if hasattr(mesh, 'nz'):
        raise Exception, "3D meshes not yet implemented"
    elif hasattr(mesh, 'ny'):
        dx = (mesh.dy, mesh.dx)
        shape = (mesh.ny, mesh.nx)
    elif hasattr(mesh, 'nx'):
        dx = (mesh.dx,)
        shape = mesh.shape
    else:
        raise Exception, "Non grid meshes can not be used for solving the FMM."
    
    phi = numerix.reshape(phi, shape)

    if LSM_SOLVER == 'lsmlib':
        from pylsmlib import computeDistanceFunction as distance
    elif LSM_SOLVER == 'skfmm':
        from skfmm import distance
    else:
        raise Exception, "Neither `lsmlib` nor `skfmm` can be found on the $PATH"

    return distance(numerix.reshape(phi, shape), dx=dx, order=order).flatten()
