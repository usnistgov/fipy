#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cahnHilliardTanh1D.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 6/15/04 {11:14:08 AM}
 # Stolen from:
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

"""

This example solves the Cahn-Hilliard equation given by:

.. raw:: latex

    $$ \\frac{\\partial \\phi}{\\partial t} = \\nabla \\cdot D \\nabla \\left( \\frac{\\partial f}{\\partial \\phi} - \\epsilon^2 \\nabla^2 \\phi \\right) $$

where the free energy functional is given by,

.. raw:: latex

    $$ f = \\frac{a^2}{2} \\phi^2 (1 - \\phi)^2 $$

The solution to the 1D problem over an infinite domain is given by,

.. raw:: latex

    $$ \\phi(x) = \\frac{1}{1 + \exp{\\left(-\\frac{a}{\\epsilon} x \\right)}} $$

Evolve the solution to equilibrium,

   >>> for step in range(steps):
   ...     it.timestep(dt = 10)

Calculate the answer,

   >>> a = Numeric.sqrt(parameters['asq'])
   >>> answer = 1 / (1 + Numeric.exp(-a * (mesh.getCellCenters()[:,0] - L / 2) / parameters['epsilon']))

Compare with the numerical solution,

   >>> Numeric.allclose(Numeric.array(var), answer, atol = 1e-2)
   1
   >>> Numeric.allclose(Numeric.array(var), answer, atol = 1e-3)
   1
   >>> Numeric.allclose(Numeric.array(var), answer, atol = 1e-4)
   0
   >>> Numeric.allclose(Numeric.array(var), answer, atol = 1e-5)
   0
"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
from fipy.iterators.iterator import Iterator
from fipy.meshes.grid2D import Grid2D
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.solvers.linearLUSolver import LinearLUSolver
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.models.cahnHilliard.cahnHilliardEquation import CahnHilliardEquation
from fipy.equations.nthOrderDiffusionEquation import NthOrderDiffusionEquation

L = 1.
nx = 100
ny = 100
steps = 100
dx = L / nx
dy = L / ny

parameters={
    'asq' : 1.0,
    'epsilon' : 1,
    'diffusionCoeff' : 1}
         

mesh = Grid2D(dx, dy, nx, ny)

##a = Numeric.sqrt(parameters['asq'])
##answer = 1 / (1 + Numeric.exp(a * (mesh.getCellCenters()[:,0] - L / 2) / parameters['epsilon']))
import random
var = CellVariable(
    name = "phase field",
    mesh = mesh,
    value = [random.random() for x in range(nx * ny)])

##var.setValue(1, cells = mesh.getCells(lambda cell: cell.getCenter()[0] > L / 2))

eqch= CahnHilliardEquation(
    var,
    parameters = parameters,
    solver = LinearPCGSolver(tolerance = 1e-15,steps = 100),
    boundaryConditions=(
    FixedFlux(mesh.getFacesRight(), 0),
    FixedFlux(mesh.getFacesLeft(), 0),
    NthOrderBoundaryCondition(mesh.getFacesLeft(), 0, 2),
    NthOrderBoundaryCondition(mesh.getFacesRight(), 0, 2),
    NthOrderBoundaryCondition(mesh.getFacesTop(), 0, 2),
    NthOrderBoundaryCondition(mesh.getFacesBottom(), 0, 2)))

it = Iterator((eqch,))

if __name__ == '__main__':
    viewer = Grid2DGistViewer(var, minVal=0., maxVal=1.0, palette = 'rainbow.gp')
    viewer.plot()
    dexp=-5

    a = Numeric.sqrt(parameters['asq'])
    answer = 1 / (1 + Numeric.exp(-a * (mesh.getCellCenters()[:,0]) / parameters['epsilon']))
    
    for step in range(steps):      
        dt=Numeric.exp(dexp)
        dt = min(10,dt)
        dexp=dexp+.5
##        print 'in loop ',step,dt
        print 'step:',step,dt
        it.timestep(dt = dt)
        diff = abs(answer - Numeric.array(var))
        maxarg = Numeric.argmax(diff)
        print 'maximum error:',diff[maxarg]
        print 'element id:',maxarg
        print 'value at element ',maxarg,' is ',var[maxarg]
        
        
        viewer.plot()
        
        
##    a = Numeric.sqrt(parameters['asq'])
##    answer = 1 / (1 + Numeric.exp(a * (mesh.getCellCenters()[:,0] - L / 2) / parameters['epsilon']))

##    from fipy.viewers.gist1DViewer import Gist1DViewer
##    view1Dvar = Gist1DViewer(vars = (var,), minVal = 0., maxVal = 1.)
##    view1Ddiff = Gist1DViewer(vars = (var - answer,), minVal = -1e-3, maxVal = 1e-3)

##    view1Dvar.plot()
##    view1Ddiff.plot()
    
##    diff = answer - Numeric.array(var)
##    maxarg = Numeric.argmax(diff)
##    print diff[maxarg]
##    print maxarg
##    aa = Numeric.argsort(diff)
##    print diff[aa[nx-1]],diff[aa[nx-2]],aa[nx-1],aa[nx-2]
##    print diff[800:840]
##    print answer[800:840]

    
    raw_input('finished')


