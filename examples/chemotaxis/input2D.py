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
    >>> accuracy = 1e-2
    >>> Numeric.allclose(KMVar, params['KM'], atol = accuracy)
    1
    >>> Numeric.allclose(TMVar, params['TM'], atol = accuracy)
    1
    >>> Numeric.allclose(TCVar, params['TC'], atol = accuracy)
    1
    >>> Numeric.allclose(P2Var, params['P2'], atol = accuracy)
    1
    >>> Numeric.allclose(P3Var, params['P3'], atol = accuracy)
    1
    >>> Numeric.allclose(KCVar, params['KC'], atol = accuracy)
    1

"""

from parameters import parameters
from fipy.meshes.grid2D import Grid2D
from fipy.variables.cellVariable import CellVariable
from fipy.equations.sourceEquation import SourceEquation
from fipy.equations.diffusionEquationWithSource import DiffusionEquationWithSource
from fipy.iterators.iterator import Iterator
from fipy.variables.variable import Variable
from fipy.boundaryConditions.fixedValue import FixedValue

params = parameters['case 2']

nx = 50
ny = 50
dx = 1.
L = nx * dx

mesh = Grid2D(nx = nx, ny = ny, dx = dx, dy = 1.)

shift = 1.

KMVar = CellVariable(mesh = mesh, value = params['KM'] * shift, hasOld = 1)
KCVar = CellVariable(mesh = mesh, value = params['KC'] * shift, hasOld = 1)
TMVar = CellVariable(mesh = mesh, value = params['TM'] * shift, hasOld = 1)
TCVar = CellVariable(mesh = mesh, value = params['TC'] * shift, hasOld = 1)
P3Var = CellVariable(mesh = mesh, value = params['P3'] * shift, hasOld = 1)
P2Var = CellVariable(mesh = mesh, value = params['P2'] * shift, hasOld = 1)
RVar = CellVariable(mesh = mesh, value = params['R'], hasOld = 1)

PN = P3Var + P2Var

KMEq = SourceEquation(KMVar,
                      scCoeff = params['chiK'] * (RVar + 1) * (1 - KCVar - KMVar.getCellVolumeAverage()),
                      spCoeff = params['lambdaK'] / (1 + PN / params['kappaK']))
TMEq = SourceEquation(TMVar,
                      scCoeff = params['chiT'] * (1 - TCVar - TMVar.getCellVolumeAverage()),
                      spCoeff = params['lambdaT'] * (KMVar + params['zetaT']))

TCEq = SourceEquation(TCVar,
                      scCoeff = params['lambdaT'] * (TMVar * KMVar).getCellVolumeAverage(),
                      spCoeff = params['lambdaTstar'])

PIP2PITP = PN / (PN / params['kappam'] + PN.getCellVolumeAverage() / params['kappac'] + 1) + params['zetaPITP']

P3Eq = DiffusionEquationWithSource(P3Var,
                                   diffusionCoeff = params['diffusionCoeff'],
                                   scCoeff = params['chi3'] * KMVar * (PIP2PITP / (1 + KMVar / params['kappa3']) + params['zeta3PITP']) + params['zeta3'],
                                   spCoeff = params['lambda3'] * (TMVar + params['zeta3T']))

P2Eq = DiffusionEquationWithSource(P2Var,
                                   diffusionCoeff = params['diffusionCoeff'],
                                   scCoeff = params['chi2'] + params['lambda3'] * params['zeta3T'] * P3Var,
                                   spCoeff = params['lambda2'] * (TMVar + params['zeta2T']))

KCEq = SourceEquation(KCVar,
                      scCoeff = params['alphaKstar'] * params['lambdaK'] * (KMVar / (1 + PN / params['kappaK'])).getCellVolumeAverage(),
                      spCoeff = params['lambdaKstar'] / (params['kappaKstar'] + KCVar))


it = Iterator((KMEq, TMEq, TCEq, P3Eq, P2Eq, KCEq))

from fipy.viewers.gnuplotViewer import GnuplotViewer
import Numeric
PNArray = Numeric.reshape(PN / PN.getCellVolumeAverage(), (50, 50))
PNViewer = GnuplotViewer(PNArray, maxVal = 2., minVal = 0., title = 'Normalised PN distribution')

KMArray = Numeric.reshape(KMVar / KMVar.getCellVolumeAverage(), (50, 50))
KMViewer = GnuplotViewer(KMArray, maxVal = 2., minVal = 0., title = 'Normalised KM distribution')

TMArray = Numeric.reshape(TMVar / TMVar.getCellVolumeAverage(), (50, 50))
TMViewer = GnuplotViewer(TMArray, maxVal = 2., minVal = 0., title = 'Normalised TM distribution')

if __name__ == '__main__':

    for i in range(100):
        print i
        it.timestep(dt = 1.)

    import Numeric

    x = mesh.getCellCenters()[:,0]
    y = mesh.getCellCenters()[:,1]

    RVar[:] = L / Numeric.sqrt((x - L / 2)**2 + (y - 2 * L)**2)
    
    for i in range(100):
        it.timestep(dt = 0.1)

    PNArray[:] = Numeric.reshape(PN / PN.getCellVolumeAverage(), (50, 50))
    PNViewer.plot('PN.ps')

    KMArray[:] = Numeric.reshape(KMVar / KMVar.getCellVolumeAverage(), (50, 50))
    KMViewer.plot('KM.ps')

    TMArray[:] = Numeric.reshape(TMVar / TMVar.getCellVolumeAverage(), (50, 50))
    TMViewer.plot('TM.ps')

    raw_input("finished")

    


