#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input1Dpoisson.py"
 #                                    created: 1/15/04 {3:45:27 PM} 
 #                                last update: 1/26/04 {6:25:33 PM} 
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
 #  2004-01-15 JEG 1.0 original
 # ###################################################################
 ##

from fivol.meshes.grid2D import Grid2D

from fivol.boundaryConditions.fixedValue import FixedValue
from fivol.boundaryConditions.fixedFlux import FixedFlux
from fivol.iterators.iterator import Iterator
from fivol.solvers.linearLUSolver import LinearLUSolver
from fivol.variables.cellVariable import CellVariable

import fivol.tools.dimensions.physicalField as physicalField

from componentVariable import ComponentVariable
from solventVariable import SolventVariable
from poissonEquation import PoissonEquation

nx = 200
dx = 0.01
L = nx * dx

parameters = {
    'mesh': {
	'nx': nx,
	'ny': 1,
	'dx': dx,
	'dy': dx
    }
}

mesh = Grid2D(
    dx = dx,
    dy = dx,
    nx = nx,
    ny = 1)

fields = {}
fields['all'] = []

fields['potential'] = CellVariable(
    mesh = mesh,
    name = 'psi',
    value = 0
    )
fields['all'] += [fields['potential']]
    
fields['substitutionals'] = []

fields['solvent'] = SolventVariable(
    mesh = mesh,
    parameters = {
	'standard potential': 1,
	'barrier height': 0
    },
    substitutionals = fields['substitutionals']
    )
	
fields['interstitials'] = (ComponentVariable(
    mesh = mesh,
    parameters = {
	'name': "e-",
	'valence': -1,
	'standard potential': 1,
	'barrier height': 0.
    }
    ),
)
fields['all'] += [fields['interstitials'][0]]

## physicalField.AddConstant(name = 'Rgas', constant = 'Nav*kB')
## physicalField.AddConstant(name = 'Faraday', constant = 'Nav*e')
physicalField.AddConstant(name = 'Rgas', constant = 1)
physicalField.AddConstant(name = 'Faraday', constant = 1)
physicalField.AddConstant(name = 'LENGTH', constant = 1)
physicalField.AddConstant(name = 'ENERGY', constant = 1)
physicalField.AddConstant(name = 'MOLARVOLUME', constant = 1)
    
poisson = PoissonEquation(
    potential = fields['potential'],
    parameters = {
	'dielectric': 1.
    },
    fields = fields,
    solver = LinearLUSolver(),
    boundaryConditions=(
	FixedValue(faces = mesh.getFacesLeft(),value = 0 * fields['potential'][0]),
	FixedFlux(faces = mesh.getFacesRight(),value = 0),
	FixedFlux(faces = mesh.getFacesTop(),value = 0),
	FixedFlux(faces = mesh.getFacesBottom(),value = 0)
    )
)

it = Iterator(equations = (poisson,))

fields['all'] += (fields['charge'],)
