#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input2D.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 9/15/05 {7:03:06 PM}
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

r"""

This example solves the Cahn-Hilliard equation given by:

.. raw:: latex

    $$ \frac{\partial \phi}{\partial t} = \nabla \cdot D \nabla
    \left( \frac{\partial f}{\partial \phi} - \epsilon^2
    \nabla^2 \phi \right) $$

where the free energy functional is given by,

.. raw:: latex

    $$ f = \frac{a^2}{2} \phi^2 (1 - \phi)^2 $$

The equation is transformed into
the following form,

.. raw:: latex

    $$ \frac{\partial \phi}{\partial t} 
    = \nabla \cdot D \frac{\partial^2 f}{\partial \phi^2} \nabla \phi - \nabla \cdot D \nabla \epsilon^2 \nabla^2 \  phi $$

This form of the equation allows the `CahnHilliardEquation` to be
constructed from a transient term, a diffusion term, and a fourth
order diffusion term. Notice that the diffusion coefficient for the
diffusion term does not always remain positive since,

.. raw:: latex

    $$ \frac{\partial^2 f}{\partial \phi^2} = a^2 (1 - 6 \phi (1 - \phi)) $$

can be less than zero and thus unstable. The fourth order diffusion
term acts to stabilize the problem.

"""
__docformat__ = 'restructuredtext'

from fipy.tools.parser import parse

numberOfElements = parse('--numberOfElements', action = 'store', type = 'int', default = 400)
numberOfSteps = parse('--numberOfSteps', action = 'store', type = 'int', default = 10)

import fipy.tools.numerix as numerix
nx = int(numerix.sqrt(numberOfElements))
ny = int(numerix.sqrt(numberOfElements))

steps = numberOfSteps

dx = 2.
dy = 2.

L = dx * nx

asq = 1.0
epsilon = 1
diffusionCoeff = 1

from fipy.meshes.grid2D import Grid2D
mesh = Grid2D(dx, dy, nx, ny)

from fipy.variables.cellVariable import CellVariable
from fipy.tools.numerix import random

var = CellVariable(name = "phase field",
                   mesh = mesh,
                   value = random.random(nx * ny))

faceVar = var.getArithmeticFaceValue()
doubleWellDerivative = asq * ( 1 - 6 * faceVar * (1 - faceVar))

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.transientTerm import TransientTerm
diffTerm2 = ImplicitDiffusionTerm(coeff = (diffusionCoeff * doubleWellDerivative,))
diffTerm4 = ImplicitDiffusionTerm(coeff = (diffusionCoeff, -epsilon**2))
eqch = TransientTerm() - diffTerm2 - diffTerm4

from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.solvers.linearLUSolver import LinearLUSolver
##solver = LinearLUSolver(tolerance = 1e-15,steps = 1000)
solver = LinearPCGSolver(tolerance = 1e-15,steps = 1000)

from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
BCs = (FixedFlux(mesh.getFacesRight(), 0),
       FixedFlux(mesh.getFacesLeft(), 0),
       NthOrderBoundaryCondition(mesh.getFacesLeft(), 0, 3),
       NthOrderBoundaryCondition(mesh.getFacesRight(), 0, 3),
       NthOrderBoundaryCondition(mesh.getFacesTop(), 0, 3),
       NthOrderBoundaryCondition(mesh.getFacesBottom(), 0, 3))

if __name__ == '__main__':

    import fipy.viewers
    viewer = fipy.viewers.make(vars = var, limits = {'datamin': 0., 'datamax': 1.0})
    viewer.plot()
    
dexp=-5

for step in range(steps):
    dt = numerix.exp(dexp)
    dt = min(100, dt)
    dexp += 0.01
    var.updateOld()
    eqch.solve(var, boundaryConditions = BCs, solver = solver, dt = dt)

    if __name__ == '__main__':
        viewer.plot()
        print 'step',step,'dt',dt
	
def _run():
    pass
            
