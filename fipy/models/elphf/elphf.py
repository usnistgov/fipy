#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "elphf.py"
 #                                    created: 12/12/03 {10:41:56 PM} 
 #                                last update: 12/8/04 {5:27:51 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from fipy.variables.cellVariable import CellVariable
from fipy.tools.dimensions import physicalField

from phaseVariable import PhaseVariable
from componentVariable import ComponentVariable
from substitutionalVariable import SubstitutionalVariable
from solventVariable import SolventVariable
 
from phaseEquation import factory as phaseFactory
from poissonEquation import factory as poissonFactory
from substitutionalEquation import factory as substitutionalFactory
from interstitialEquation import factory as interstitialFactory

constant = {}

def addScales(mesh, parameters):
    
    def getScale(unit, parameter = None, value = None):
	if parameter is not None and parameters.has_key(parameter):
	    value = parameters[parameter]
	elif value is None:
	    value = 1
	value = physicalField.PhysicalField(value)
	if not value.getUnit().isDimensionless():
	    value.convertToUnit(unit)
	return value
	
    constant['TIME'] = getScale(unit = 's', parameter = 'time step duration')
    constant['LENGTH'] = getScale(unit = 'm', value = mesh.getScale())
    mesh.setScale(1)
    
    constant['MOLARVOLUME'] = getScale(unit = 'm**3/mol', parameter = 'substitutional molar volume', value = '1 m**3/mol')
    TEMPERATURE = getScale(unit = 'K', parameter = 'temperature', value = '298 K')
    Rgas = physicalField.PhysicalField(value = '1 Nav*kB')
    constant['ENERGY'] = getScale(unit = 'J/mol', value = Rgas * TEMPERATURE)
    constant['Faraday'] = physicalField.PhysicalField(value = '1 Nav*e')
    constant['POTENTIAL'] = getScale(unit = 'V', value = constant['ENERGY'] / constant['Faraday'])
    
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
	    'permittivity': 1.
	}
	addPotentialToAll = 0
	
    if not parameters['potential'].has_key('permittivity'):
	parameters['potential']['permittivity'] = parameters['potential']['dielectric'] * physicalField.PhysicalField(1, "eps0")

    parameters['potential']['var'] = CellVariable(
	mesh = mesh,
	name = parameters['potential']['name'],
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
		value = value
	    )
	    
	    component['var'].parameters = component
	    fields['interstitials'] += (component['var'],)
	    fields['all'] += [component['var']]
    
    fields['substitutionals'] = ()
    
    if not parameters.has_key('solvent'):
	parameters['solvent'] = {
	    'standard potential': 0,
	    'barrier height': 0
	}
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
		solventParameters = parameters['solvent'],
		value = value
	    )
	    
	    component['var'].parameters = component
	    fields['substitutionals'] += (component['var'],)
	    fields['all'] += [component['var']]
	    
	addSolventToAll = 1

    parameters['solvent']['var'] = SolventVariable(
	mesh = mesh,
	parameters = parameters['solvent'],
	substitutionals = fields['substitutionals']
	)
	
    fields['solvent'] = parameters['solvent']['var']
    fields['solvent'].parameters = parameters['solvent']
    
    if addSolventToAll:
	fields['all'] += [fields['solvent']]
	
    fields['charge'] = fields['solvent'].getValence()
    for component in list(fields['interstitials']) + list(fields['substitutionals']):
	fields['charge'] = fields['charge'] + component * component.getValence() #.getOld()

    return fields

def makeEquations(fields, parameters):
    fields['potential'].equation = poissonFactory.make(fields, parameters['potential'])
    fields['phase'].equation = phaseFactory.make(fields, parameters)
    
    for component in fields['substitutionals']:
	component.equation = substitutionalFactory.make(component, fields)

    for component in fields['interstitials']:
	component.equation = interstitialFactory.make(component, fields)

