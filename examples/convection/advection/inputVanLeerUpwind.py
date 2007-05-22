#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "inputExplicitUpwind.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 5/18/06 {8:40:11 PM} 
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

r""" 

This example demonstrates the use of the `VanLeerConvectionTerm` as
defined by http://www.gre.ac.uk/~physica/phy2.12/theory/node173.htm

In this example a square wave is advected. The Van Leer discretization
should in theory do a good job of preserving the shape of the
wave. This may or may not be happening in this case. This example
needs further testing

The test case is mainly to check that the periodic mesh is working
correctly. We advect the wave on different meshes one periodic and one
non-periodic but twice as long. The results are then compared. The
periodic wave wraps around the mesh.

    >>> newVar2 = var2.copy()

    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> for step in range(steps):
    ...	    eq1.solve(var = var1, dt = dt, solver = LinearLUSolver())
    ...     eq2.solve(var = var2, dt = dt, solver = LinearLUSolver())

    
    >>> newVar2[:nx / 4] = var2[nx / 4:]
    >>> newVar2[nx / 4:] = var2[:nx / 4]
    >>> print newVar2.allclose(var1[nx / 4:3 * nx / 4], atol = 1e-6)
    1

Currently after 20 steps the wave has lost 23% of its height. Van Leer
should do better than this.
    
    >>> var1.max() > 0.77
    1
    
"""

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
     
L = 20.
nx = 40
dx = L / nx
cfl = 0.5
velocity = 1.0
dt = cfl * dx / velocity

steps = int(L /  4. / dt / velocity)

from fipy.meshes.grid1D import Grid1D
mesh = Grid1D(dx = dx, nx = nx)

from fipy.meshes.periodicGrid1D import PeriodicGrid1D
periodicMesh = PeriodicGrid1D(dx = dx, nx = nx / 2)

startingArray = numerix.zeros(nx, 'd')
startingArray[2 * nx / 10: 3 * nx / 10] = 1. 

from fipy.variables.cellVariable import CellVariable
var1 = CellVariable(
    name = "non-periodic",
    mesh = mesh,
    value = startingArray)

var2 = CellVariable(
    name = "periodic",
    mesh = periodicMesh,
    value = startingArray[:nx / 2])

from fipy.terms.transientTerm import TransientTerm
from fipy.terms.vanLeerConvectionTerm import VanLeerConvectionTerm
eq1 = TransientTerm() - VanLeerConvectionTerm(coeff = (-velocity,))
eq2 = TransientTerm() - VanLeerConvectionTerm(coeff = (-velocity,))

if __name__ == '__main__':

    import fipy.viewers
    viewer1 = fipy.viewers.make(vars=var1)
    viewer2 = fipy.viewers.make(vars=var2)
    viewer1.plot()
    viewer2.plot()
    from fipy.solvers.linearLUSolver import LinearLUSolver

    newVar2 = var2.copy()

    for step in range(steps):
        eq1.solve(var = var1, dt = dt, solver = LinearLUSolver())
        eq2.solve(var = var2, dt = dt, solver = LinearLUSolver())
        viewer1.plot()
        viewer2.plot()

    newVar2[:nx / 4] = var2[nx / 4:]
    newVar2[nx / 4:] = var2[:nx / 4]

    print 'maximum absolute difference between periodic and non-periodic grids:',numerix.max(abs(var1[nx / 4:3 * nx / 4] - newVar2))

    raw_input('finished')
