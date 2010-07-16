#!/usr/bin/env python

## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "collectedDiffusionTerm.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # #############################################################################
 ##

from fipy.terms.diffusionTerm import DiffusionTerm

class _CollectedDiffusionTerm(DiffusionTerm):
    def __init__(self):
        self.orders = [None, None]

    def __iadd__(self, other):
        if not isinstance(other, DiffusionTerm):
            return DiffusionTerm.__add__(self, other)
        elif isinstance(other, _CollectedDiffusionTerm):
            for term in other.orders:
                self += term
        elif other.order == 0:
            if self.orders[0] is None:
                self.orders[0] = other
        elif other.order == 2:
            # index is order / 2
            if self.orders[1] is None:
                self.orders[1] = other
            else:
                self.orders[1] += other
        else:
            self.orders.append(other)
            
        return self
                
    def __add__(self, other):
        dup = self.copy()
        dup += other
        
        return dup
        
    def _buildMatrix(self, var, SparseMatrix, boundaryConditions, dt, equation=None):
        from fipy.tools import numerix
        from fipy.matrices.sparseMatrix import _SparseMatrix
        
        N = len(var)
        RHSvector = numerix.zeros((N,),'d')
        matrix = SparseMatrix(mesh=var.getMesh())

        for term in self.orders:
            if term is not None:
                termMatrix, termRHSvector = term._buildMatrix(var, SparseMatrix,
                                                              boundaryConditions, 
                                                              dt=dt, equation=equation)

                matrix += termMatrix
                RHSvector += termRHSvector
        
        return (matrix, RHSvector)

    def copy(self):
        dup = _CollectedDiffusionTerm()
        dup.orders = self.orders
        return dup
        
    def _getDefaultSolver(self, solver, *args, **kwargs):
        for term in self.orders:
            if term is not None:
                solver = term._getDefaultSolver(solver, *args, **kwargs)
                if solver is not None:
                    return solver
                
        return None
        
    def __repr__(self):
        reprs = []
        for term in self.orders:
            if term is not None:
                reprs.append(repr(term))
                
        return " + ".join(reprs)

    def __neg__(self):
        r"""
         Negate a `Term`.

           >>> -(DiffusionTerm(coeff=1.) - DiffusionTerm(coeff=(2., -3.)))
           DiffusionTerm(coeff=[-1.0]) + DiffusionTerm(coeff=[2.0, -3.0])

        """
        dup = _CollectedDiffusionTerm()
        
        dup.orders = []
        for term in self.orders:
            if term is None:
                dup.orders.append(term)
            else:
                dup.orders.append(-term)
                
        return dup
        
    def _getGeomCoeff(self, mesh):
        # index is order / 2
        if self.orders[1] is not None:
            return self.orders[1]._getGeomCoeff(mesh)
        else:
            return None


    def _test(self):
        """
        Introduced for __iadd__ bug.

           >>> from fipy import *
           >>> TransientTerm() == DiffusionTerm((0, 0)) + DiffusionTerm() + ImplicitSourceTerm()
           TransientTerm(coeff=1.0) + DiffusionTerm(coeff=[-1.0]) + DiffusionTerm(coeff=[0, 0]) + ImplicitSourceTerm(coeff=-(0.0)) == 0
           
        """
        pass
    
def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()

        
