#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:02:23 PM} 
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

This example solves the steady-state convection-diffusion equation as described in
`./examples/diffusion/convection/exponential1D/inpuy.py` but uses the
`PowerLawConvectionTerm` rather than the `ExponentialConvectionTerm` instatiator.

NOTE: This problem does not appear to converge with a triangular mesh (last tested 7/8/2004). It is believed that this may be bue to the PCGSolver.

The analytical solution test for this problem is given by:

   >>> axis = 0
   >>> x = mesh.getCellCenters()[:,axis]
   >>> AA = -sourceCoeff * x / convCoeff[axis]
   >>> BB = 1. + sourceCoeff * L / convCoeff[axis]
   >>> import Numeric
   >>> CC = 1. - Numeric.exp(-convCoeff[axis] * x / diffCoeff)
   >>> DD = 1. - Numeric.exp(-convCoeff[axis] * L / diffCoeff)
   >>> analyticalArray = AA + BB * CC / DD
   >>> Numeric.allclose(analyticalArray, Numeric.array(var), rtol = 1e-2, atol = 1e-2) 
   1
   
"""
     
from fipy.meshes.numMesh.tri2D import Tri2D
from fipy.equations.stdyConvDiffScEquation import SteadyConvectionDiffusionScEquation
from fipy.solvers.linearCGSSolver import LinearCGSSolver
from fipy.iterators.iterator import Iterator
from fipy.variables.cellVariable import CellVariable
from fipy.terms.scSourceTerm import ScSourceTerm
from fipy.viewers.pyxviewer import PyxViewer
from fipy.terms.powerLawConvectionTerm import PowerLawConvectionTerm
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux

valueLeft = 0.
valueRight = 1.
L = 10.
nx = 1000
ny = 1
diffCoeff = 1.
convCoeff = (10.,0.)
sourceCoeff = 0.

mesh = Tri2D(L / nx, L / ny, nx, ny)

var = CellVariable(
    name = "solution variable",
    mesh = mesh,
    value = valueLeft)

boundaryConditions = (
    FixedValue(mesh.getFacesLeft(), valueLeft),
    FixedValue(mesh.getFacesRight(), valueRight),
    FixedFlux(mesh.getFacesTop(), 0.),
    FixedFlux(mesh.getFacesBottom(), 0.)
    )

        
eq = SteadyConvectionDiffusionScEquation(
    var = var,
    diffusionCoeff = diffCoeff,
    convectionCoeff = convCoeff,
    sourceCoeff = sourceCoeff,
    solver = LinearCGSSolver(tolerance = 1.e-15, steps = 2000),
    convectionScheme = PowerLawConvectionTerm,
    boundaryConditions = boundaryConditions
    )

it = Iterator((eq,))

it.timestep()

if __name__ == '__main__':
    viewer = PyxViewer(var)
    print var
    viewer.plot()
    raw_input('finished')
