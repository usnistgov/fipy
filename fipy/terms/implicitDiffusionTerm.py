#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "implicitDiffusionTerm.py"
 #                                    created: 11/28/03 {10:07:06 AM} 
 #                                last update: 12/6/04 {4:54:19 PM} 
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


from fipy.terms.diffusionTerm import DiffusionTerm

class ImplicitDiffusionTerm(DiffusionTerm):
    """
	>>> from fipy.meshes.grid2D import Grid2D
	>>> mesh = Grid2D(nx = 2)
	
	>>> from fipy.variables.cellVariable import CellVariable
	>>> var = CellVariable(mesh = mesh)
	
	>>> term = ImplicitDiffusionTerm()
	
	>>> from fipy.boundaryConditions.fixedValue import FixedValue
	>>> bcs = (FixedValue(faces = mesh.getFacesLeft(), value = 0), 
	...     FixedValue(faces = mesh.getFacesRight(), value = 1))
	>>> term.solve(var = var, boundaryConditions = bcs)
	>>> print var
	[ 0.25, 0.75,]
	
	>>> eq = term + ImplicitDiffusionTerm()
	>>> bcs = (FixedValue(faces = mesh.getFacesLeft(), value = 3), 
	...     FixedValue(faces = mesh.getFacesRight(), value = 2))
	>>> eq.solve(var = var, boundaryConditions = bcs)
	>>> print var
	[ 2.75, 2.25,]
    """
    def getWeight(self):
	return {
	    'implicit':{
		'cell 1 diag':     1, 
		'cell 1 offdiag': -1, 
		'cell 2 diag':     1, 
		'cell 2 offdiag': -1
	    }
	}

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

