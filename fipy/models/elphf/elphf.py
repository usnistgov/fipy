#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "elphf.py"
 #                                    created: 12/12/03 {10:41:56 PM} 
 #                                last update: 12/29/03 {1:22:30 PM} 
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

 
from variables.cellVariable import CellVariable
from phaseVariable import PhaseVariable
from componentVariable import ComponentVariable
from substitutionalVariable import SubstitutionalVariable
from solventVariable import SolventVariable
 
from phaseEquation import PhaseEquation
from poissonEquation import PoissonEquation
from substitutionalEquation import SubstitutionalEquation
from interstitialEquation import InterstitialEquation

from solvers.linearCGSSolver import LinearCGSSolver
from solvers.linearLUSolver import LinearLUSolver
from solvers.linearGMRESSolver import LinearGMRESSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator

def makeFields(mesh, parameters):
    fields = {}
    
    fields['all'] = []
    
    parameters['phase']['var'] = PhaseVariable(
	name = parameters['phase']['name'],
	mesh = mesh,
	value = parameters['phase']['initial'][0],
	)
	
    fields['phase'] = parameters['phase']['var']
    
    fields['all'] += [fields['phase']]
    
    parameterList = [parameters['phase']]
    
    if not parameters.has_key('potential'):
	parameters['potential'] = {
	    'name': "psi",
	    'permittivity': 1.
	}

    parameters['potential']['var'] = CellVariable(
	mesh = mesh,
	name = 'psi',
	)
    fields['potential'] = parameters['potential']['var']
	
    fields['all'] += [fields['potential']]
    
    fields['interstitials'] = ()
    
    if parameters.has_key('interstitials'):
	for component in parameters['interstitials']:
	    component['var'] = ComponentVariable(
		mesh = mesh,
		parameters = component,
		value = component['initial'][0],
		)
	    
	    fields['interstitials'] += (component['var'],)
	    fields['all'] += [component['var']]
	    
	    parameterList += list(parameters['interstitials'])
    
    fields['substitutionals'] = ()
    
    if parameters.has_key('substitutionals'):
	for component in parameters['substitutionals']:
	    component['var'] = SubstitutionalVariable(
		mesh = mesh,
		parameters = component,
		solventParameters = parameters['solvent'],
		value = component['initial'][0],
		)
	    
	    fields['substitutionals'] += (component['var'],)
	    fields['all'] += [component['var']]

	    parameterList += list(parameters['substitutionals'])
	    
    parameters['solvent']['var'] = SolventVariable(
	mesh = mesh,
	parameters = parameters['solvent'],
	substitutionals = fields['substitutionals']
	)
	
    fields['solvent'] = parameters['solvent']['var']
    fields['all'] += [fields['solvent']]
	
    # set initial conditions
    for field in parameterList:
	for init in field['initial'][1:]:
	    setCells = mesh.getCells(init['func'])
	    field['var'].setValue(init['value'],setCells)
    
    return fields


def makeIterator(mesh, fields, parameters, maxSweeps = 1):
    equations = (PhaseEquation(
	phase = fields['phase'],
	timeStepDuration = parameters['time step duration'],
	fields = fields,
	phaseMobility = parameters['phase']['mobility'],
	phaseGradientEnergy = parameters['phase']['gradient energy'],
	solver = LinearLUSolver(),
	boundaryConditions=(
# 	    FixedValue(faces = mesh.getFacesLeft(),value = 1.),
# 	    FixedValue(faces = mesh.getFacesRight(),value = 0.),
	    FixedFlux(faces = mesh.getFacesLeft(),value = 0.),
	    FixedFlux(faces = mesh.getFacesRight(),value = 0.),
	    FixedFlux(faces = mesh.getFacesTop(),value = 0.),
	    FixedFlux(faces = mesh.getFacesBottom(),value = 0.)
	)
    ),)
    
    equations += (PoissonEquation(
	    potential = fields['potential'],
	    parameters = parameters['potential'],
	    fields = fields,
	    solver = LinearLUSolver(),
	    boundaryConditions=(
# 		FixedValue(faces = mesh.getFacesLeft(),value = 1.),
# 		FixedValue(faces = mesh.getFacesRight(),value = 0.),
		FixedValue(faces = mesh.getFacesLeft(),value = 0.),
		FixedFlux(faces = mesh.getFacesRight(),value = 0.),
		FixedFlux(faces = mesh.getFacesTop(),value = 0.),
		FixedFlux(faces = mesh.getFacesBottom(),value = 0.)
	    )
	),
    )
    
    for component in fields['substitutionals']:
	eq = SubstitutionalEquation(
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
	
    for component in fields['interstitials']:
	eq = InterstitialEquation(
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


