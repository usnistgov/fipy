#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "elphf.py"
 #                                    created: 12/12/03 {10:41:56 PM} 
 #                                last update: 12/22/03 {5:01:22 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from concentrationEquation import ConcentrationEquation
from solventVariable import SolventVariable

from solvers.linearCGSSolver import LinearCGSSolver
from solvers.linearLUSolver import LinearLUSolver
from solvers.linearGMRESSolver import LinearGMRESSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator

def makeIterator(mesh, fields, parameters, maxSweeps = 1):
    equations = ()
    
    fields['solvent'] = 1.
    for component in fields['substitutionals']:
	fields['solvent'] = fields['solvent'] - component#.getOld()    
    for component in fields['substitutionals']:
	eq = ConcentrationEquation(
	    Cj = component,
	    timeStepDuration = parameters['time step duration'],
	    fields = fields,
	    diffusivity = parameters['diffusivity'],
	    solver = LinearLUSolver(),
# 	    solver = LinearGMRESSolver(
# 		 tolerance = 1.e-15, 
# 		 steps = 1000
# 	    ),
# 	   solver = LinearCGSSolver(
# 		tolerance = 1.e-15, 
# 		steps = 1000
# 	   ),
	    boundaryConditions=(
# 		FixedValue(faces = mesh.getFacesLeft(),value = parameters['valueLeft']),
# 		FixedValue(faces = mesh.getFacesRight(),value = parameters['valueRight']),
		FixedFlux(faces = mesh.getFacesLeft(),value = 0.),
		FixedFlux(faces = mesh.getFacesRight(),value = 0.),
		FixedFlux(faces = mesh.getFacesTop(),value = 0.),
		FixedFlux(faces = mesh.getFacesBottom(),value = 0.)
	    )
	)
	equations += (eq,)
	
    return Iterator(equations = equations, maxSweeps = maxSweeps)


