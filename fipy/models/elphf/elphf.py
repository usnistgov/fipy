#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "elphf.py"
 #                                    created: 12/12/03 {10:41:56 PM} 
 #                                last update: 1/14/04 {2:37:20 PM} 
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

from tools.dimensions import physicalField

def addScales(mesh, parameters):
    
    def addScale(scale, unit, parameter = None, value = None):
	if parameters.has_key('scale') and parameters['scale'].has_key(scale):
	    value = parameters['scale'][scale]
	else:
	    if parameter is not None and parameters.has_key(parameter):
		value = parameters[parameter]
	    elif value is None:
		value = 1
	physicalField.AddConstant(scale, physicalField.NonDimOrUnits(value, unit))
	
    addScale(scale = 'TIME', unit = 's', parameter = 'time step duration')
    addScale(scale = 'LENGTH', unit = 'm', value = mesh.getScale())
    mesh.setScale(1)
    addScale(scale = 'MOLARVOLUME', unit = 'm**3/mol', parameter = 'substitutional molar volume', value = '1 m**3/mol')
    addScale(scale = 'TEMPERATURE', unit = 'K', parameter = 'temperature', value = '298 K')
    physicalField.AddConstant('Rgas', 'Nav*kB')
    addScale(scale = 'ENERGY', unit = 'J/mol', value = "Rgas * TEMPERATURE")
    physicalField.AddConstant('Faraday', 'Nav*e')
    addScale(scale = 'POTENTIAL', unit = 'V', value = "ENERGY / Faraday")
    
def makeFields(mesh, parameters):
    addScales(mesh, parameters)

    fields = {}
    
    fields['all'] = []
    
    addPhaseToAll = 1
    if not parameters.has_key('phase'):
	parameters['phase'] = {
	    'name': "xi",
	    'mobility': 1,
	    'gradient energy': 1,
	    'value': 1
	}
	addPhaseToAll = 0
	
    parameters['phase']['var'] = PhaseVariable(
	name = parameters['phase']['name'],
	mesh = mesh,
## 	value = physicalField.PhysicalField(parameters['phase']['value']),
 	value = physicalField.Scale(parameters['phase']['value'], 1),
	)
	
    fields['phase'] = parameters['phase']['var']
    fields['phase'].parameters = parameters['phase']
    
    if addPhaseToAll:
	fields['all'] += [fields['phase']]
    
    addPotentialToAll = 1
    if not parameters.has_key('potential'):
	parameters['potential'] = {
	    'name': "psi",
	    'dielectric': 1.
	}
	addPotentialToAll = 0

    parameters['potential']['var'] = CellVariable(
	mesh = mesh,
	name = 'psi',
	value = 0
## 	value = "0 V"
	)
    fields['potential'] = parameters['potential']['var']
    fields['potential'].parameters = parameters['potential']
	
    if addPotentialToAll:
	fields['all'] += [fields['potential']]
    
    fields['interstitials'] = ()
    
    if parameters.has_key('interstitials'):
	for component in parameters['interstitials']:
	    if component.has_key('value'):
		value = component['value']
	    else:
		value = 0
	    component['var'] = ComponentVariable(
		mesh = mesh,
		parameters = component,
		systemParameters = parameters,
		value = value,
		)
	    
	    component['var'].parameters = component
	    fields['interstitials'] += (component['var'],)
	    fields['all'] += [component['var']]
    
    fields['substitutionals'] = ()
    
    addSolventToAll = 0
    if parameters.has_key('substitutionals'):
	for component in parameters['substitutionals']:
	    if component.has_key('value'):
		value = component['value']
	    else:
		value = 0
	    component['var'] = SubstitutionalVariable(
		mesh = mesh,
		parameters = component,
		systemParameters = parameters,
		value = value,
		)
	    
	    component['var'].parameters = component
	    fields['substitutionals'] += (component['var'],)
	    fields['all'] += [component['var']]
	addSolventToAll = 1

    parameters['solvent']['var'] = SolventVariable(
	mesh = mesh,
	parameters = parameters['solvent'],
	systemParameters = parameters,
	substitutionals = fields['substitutionals']
	)
	
    fields['solvent'] = parameters['solvent']['var']
    fields['solvent'].parameters = parameters['solvent']
    
    if addSolventToAll:
	fields['all'] += [fields['solvent']]
	
    return fields


def makeIterator(mesh, fields, parameters):
    timeStepDuration = physicalField.Scale(parameters['time step duration'], "TIME")
##     timeStepDuration = physicalField.PhysicalField(parameters['time step duration'])
    
##     gradientEnergy = physicalField.PhysicalField(parameters['phase']['gradient energy'])
##     mobility = physicalField.PhysicalField(parameters['phase']['mobility'])
##     flux = 0 * (fields['phase'][0] / mesh.getPhysicalShape()[0] ) * gradientEnergy * mobility
    flux = 0
    equations = (PhaseEquation(
	phase = fields['phase'],
	timeStepDuration = timeStepDuration,
	fields = fields,
## 	phaseMobility = mobility,
## 	phaseGradientEnergy = gradientEnergy,
 	phaseMobility = physicalField.Scale(parameters['phase']['mobility'],"MOLARVOLUME/ENERGY/TIME"),
 	phaseGradientEnergy = physicalField.Scale(parameters['phase']['gradient energy'],"LENGTH**2*ENERGY/MOLARVOLUME"),
	solver = LinearLUSolver(),
	boundaryConditions=(
# 	    FixedValue(faces = mesh.getFacesLeft(),value = 1.),
# 	    FixedValue(faces = mesh.getFacesRight(),value = 0.),
	    FixedFlux(faces = mesh.getFacesLeft(),value = flux),
	    FixedFlux(faces = mesh.getFacesRight(),value = flux),
	    FixedFlux(faces = mesh.getFacesTop(),value = flux),
	    FixedFlux(faces = mesh.getFacesBottom(),value = flux)
	)
    ),)
    
##     flux = 0 * (fields['potential'][0] / mesh.getPhysicalShape()[0] ) * "1 eps0"
    flux = 0
    equations += (PoissonEquation(
	    potential = fields['potential'],
	    parameters = parameters['potential'],
	    fields = fields,
	    solver = LinearLUSolver(),
	    boundaryConditions=(
# 		FixedValue(faces = mesh.getFacesLeft(),value = 1.),
# 		FixedValue(faces = mesh.getFacesRight(),value = 0.),
		FixedValue(faces = mesh.getFacesLeft(),value = 0 * fields['potential'][0]),
		FixedFlux(faces = mesh.getFacesRight(),value = flux), # "0 eps0*V/m"
		FixedFlux(faces = mesh.getFacesTop(),value = flux),
		FixedFlux(faces = mesh.getFacesBottom(),value = flux)
	    )
	),
    )
    
    for component in fields['substitutionals']:
## 	flux = 0 * (component[0] / mesh.getPhysicalShape()[0] ) * component.diffusivity
	flux = 0
	eq = SubstitutionalEquation(
	    Cj = component,
	    timeStepDuration = timeStepDuration,
	    fields = fields,
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
		FixedFlux(faces = mesh.getFacesLeft(),value = flux), # "0. m/s"
		FixedFlux(faces = mesh.getFacesRight(),value = flux),
		FixedFlux(faces = mesh.getFacesTop(),value = flux),
		FixedFlux(faces = mesh.getFacesBottom(),value = flux)
	    )
	)
	equations += (eq,)
	
    for component in fields['interstitials']:
## 	flux = 0 * (component[0] / mesh.getPhysicalShape()[0] ) * component.diffusivity
	flux = 0
	eq = InterstitialEquation(
	    Cj = component,
	    timeStepDuration = timeStepDuration,
	    fields = fields,
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
		FixedFlux(faces = mesh.getFacesLeft(),value = flux), # "0. m/s"
		FixedFlux(faces = mesh.getFacesRight(),value = flux),
		FixedFlux(faces = mesh.getFacesTop(),value = flux),
		FixedFlux(faces = mesh.getFacesBottom(),value = flux)
	    )
	)
	equations += (eq,)
	
    return Iterator(equations = equations)


