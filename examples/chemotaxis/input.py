#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/7/04 {8:23:02 AM} { 5:14:21 PM}
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

"""

Input file for chemotaxis modeling.

Here are some test cases for the model.

    >>> for i in range(300):
    ...     it.timestep(dt = 0.1)
    >>> import Numeric
    >>> Numeric.allclose(KMVar, params['KM'], atol = 1e-5)
    1
    >>> Numeric.allclose(TMVar, params['TM'], atol = 1e-4)
    1
    >>> Numeric.allclose(TCVar, params['TC'], atol = 1e-4)
    1

"""

from parameters import parameters
from fipy.meshes.grid2D import Grid2D
from fipy.variables.cellVariable import CellVariable
from fipy.equations.sourceEquation import SourceEquation
from fipy.iterators.iterator import Iterator

params = parameters['case 1']

mesh = Grid2D(nx = 10, ny = 1, dx = 1., dy = 1.)

KMVar = CellVariable(mesh = mesh, value = 0.2)
KCVar = CellVariable(mesh = mesh, value = params['KC'])
PNVar = CellVariable(mesh = mesh, value = params['P3'] + params['P2'])
TMVar = CellVariable(mesh = mesh, value = 0.2)
TCVar = CellVariable(mesh = mesh, value = 0.4)

KMEq = SourceEquation(KMVar,
                      scCoeff = params['chiK'] * (params['R'] + 1) * (1 - KCVar - KMVar.getCellVolumeAverage()),
                      spCoeff = params['lambdaK'] / (1 + PNVar / params['kappaK']))

TMEq = SourceEquation(TMVar,
                      scCoeff = params['chiT'] * (1 - TCVar - TMVar.getCellVolumeAverage()),
                      spCoeff = params['lambdaT'] * (KMVar + params['zetaT']))

TCEq = SourceEquation(TCVar,
                      scCoeff = params['lambdaT'] * (TMVar * KMVar).getCellVolumeAverage(),
                      spCoeff = params['lambdaTstar'])

it = Iterator((KMEq, TMEq, TCEq))

if __name__ == '__main__':

    for i in range(200):
        it.timestep(dt = 0.1)

    print KMVar, TMVar, TCVar
                      
                                                  
