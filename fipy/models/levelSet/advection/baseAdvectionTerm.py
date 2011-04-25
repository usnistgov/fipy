#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "baseAdvectionTerm.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import MA

from fipy.terms.nonDiffusionTerm import _NonDiffusionTerm

class _BaseAdvectionTerm(_NonDiffusionTerm):

    def __init__(self, coeff = None):
        if self.__class__ is _BaseAdvectionTerm:
            raise NotImplementedError, "can't instantiate abstract base class"
                                                
        _NonDiffusionTerm.__init__(self)
        self.geomCoeff = coeff

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, equation=None, transientGeomCoeff=None, diffusionGeomCoeff=None):

        oldArray = var.old

        mesh = var.mesh
        NCells = mesh.numberOfCells
        NCellFaces = mesh._maxFacesPerCell

        cellValues = numerix.repeat(oldArray[numerix.newaxis, ...], NCellFaces, axis = 0)

        cellIDs = numerix.repeat(numerix.arange(NCells)[numerix.newaxis, ...], NCellFaces, axis = 0)
        cellToCellIDs = mesh._cellToCellIDs

        if NCells > 0:
            cellToCellIDs = MA.where(MA.getmask(cellToCellIDs), cellIDs, cellToCellIDs) 

            adjacentValues = numerix.take(oldArray, cellToCellIDs)

            differences = self._getDifferences(adjacentValues, cellValues, oldArray, cellToCellIDs, mesh)
            differences = MA.filled(differences, 0)

            minsq = numerix.sqrt(numerix.sum(numerix.minimum(differences, numerix.zeros((NCellFaces, NCells)))**2, axis=0))
            maxsq = numerix.sqrt(numerix.sum(numerix.maximum(differences, numerix.zeros((NCellFaces, NCells)))**2, axis=0))

            coeff = numerix.array(self._getGeomCoeff(mesh))

            coeffXdiffereneces = coeff * ((coeff > 0.) * minsq + (coeff < 0.) * maxsq)
        else:
            coeffXdiffereneces = 0.

        return (var, SparseMatrix(mesh=var.mesh), -coeffXdiffereneces * mesh.cellVolumes)

    def _getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs, mesh):
        return (adjacentValues - cellValues) / mesh._cellToCellDistances
        
    def _getDefaultSolver(self, solver, *args, **kwargs):
        if _NonDiffusionTerm._getDefaultSolver(self, solver, *args, **kwargs) is not None:
            raise AssertionError, 'An alternate _getDefaultSolver() is defined in a base class'
        
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)

        import fipy.solvers.solver
        if fipy.solvers.solver == 'trilinos' or fipy.solvers.solver == 'no-pysparse':
            from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner
            from fipy.solvers.trilinos.linearGMRESSolver import LinearGMRESSolver
            return solver or LinearGMRESSolver(precon=JacobiPreconditioner(), *args, **kwargs)
        elif fipy.solvers.solver == 'pyamg':
            from fipy.solvers.pyAMG.linearGeneralSolver import LinearGeneralSolver
            return solver or LinearGeneralSolver(*args, **kwargs)
        else:
            from fipy.solvers import DefaultAsymmetricSolver
            return solver or DefaultAsymmetricSolver(*args, **kwargs)

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
        
        

    


        

    

    
      

        
        
        
        
