#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/13/04 {3:14:14 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
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
A simple 1D three-component diffusion problem to test the
`ConcentrationEquation` element of ElPhF. We iterate to equilibrium:

    >>> for step in range(40):
    ...     it.timestep(dt = parameters['time step duration'])
    
Verify that the concentrations have become uniform

    >>> fields['substitutionals'][0].allclose(0.45, rtol = 1e-7, atol = 1e-7)
    1
    >>> fields['substitutionals'][1].allclose(0.45, rtol = 1e-7, atol = 1e-7)
    1
"""
__docformat__ = 'restructuredtext'

## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.gist1DViewer import Gist1DViewer
from fipy.iterators.iterator import Iterator

import fipy.models.elphf.elphf as elphf

## from elphfIterator import ElPhFIterator

nx = 40
dx = 1.
L = nx * dx
mesh = Grid2D(
    dx = dx,
    dy = 1.,
    nx = nx,
    ny = 1)
    
parameters = {
    'time step duration': 10000,
    'solvent': {
	'standard potential': 0.,
	'barrier height': 0.
    }
}

parameters['substitutionals'] = (
    {
	'name': "c1",
	'diffusivity': 1.,
	'standard potential': 1.,
	'barrier height': 1.
    },
    {
	'name': "c2",
	'diffusivity': 1.,
	'standard potential': 1.,
	'barrier height': 1.
    }
)

fields = elphf.makeFields(mesh = mesh, parameters = parameters)

setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2)
fields['substitutionals'][0].setValue(0.3)
fields['substitutionals'][0].setValue(0.6,setCells)
fields['substitutionals'][1].setValue(0.6)
fields['substitutionals'][1].setValue(0.3,setCells)

equations = elphf.makeEquations(
    mesh = mesh, 
    fields = fields, 
    parameters = parameters
)

it = Iterator(equations = equations)

if __name__ == '__main__':
    viewer = Gist1DViewer(vars = fields['substitutionals'])

    viewer.plot()
	
    raw_input("press <return> to start...")

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    # it.timestep(50)
    # 
    # for timestep in range(5):
    #     it.sweep(50)
    #     it.advanceTimeStep


    for i in range(40):
	it.timestep(dt = parameters['time step duration'])
    #     raw_input()
	
	viewer.plot()
	
    # profile.stop()
	    
    raw_input("finished")

