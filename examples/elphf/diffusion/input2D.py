#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input2D.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 7/29/04 {9:30:12 AM} 
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

""" 
The same three-component diffusion problem as introduced in `input1D.py`,
but in 2D

    >>> for step in range(40):
    ...     it.timestep(dt = parameters['time step duration'])
    
Verify that the concentrations have become uniform

    >>> fields['substitutionals'][0].allclose(0.45, rtol = 1e-7, atol = 1e-7)
    1
    >>> fields['substitutionals'][1].allclose(0.45, rtol = 1e-7, atol = 1e-7)
    1
"""
__docformat__ = 'restructuredtext'
 
from fipy.tools.profiler.profiler import Profiler
from fipy.tools.profiler.profiler import calibrate_profiler

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.iterators.iterator import Iterator

import fipy.models.elphf.elphf as elphf

nx = 40
dx = 1.
L = nx * dx

## fudge = calibrate_profiler(10000)
## profile = Profiler('profile', fudge=fudge)

mesh = Grid2D(
    dx = dx,
    dy = 1.,
    nx = nx,
    ny = 40)

## profile.stop()

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
    viewers = [Grid2DGistViewer(var = field) for field in fields['all']]

    for viewer in viewers:
	viewer.plot()
	
    raw_input()

    it.timestep(1)
    
    fudge = calibrate_profiler(10000)
    profile = Profiler('profile', fudge=fudge)
    
    for i in range(10):
	it.timestep(1)
	
## 	for viewer in viewers:
## 	    viewer.plot()
	
    profile.stop()
	    
    raw_input()

