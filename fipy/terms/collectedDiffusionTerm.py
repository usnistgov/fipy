## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "collectedDiffusionTerm.py"
 #                                     created: 3/28/07 {8:44:05 AM}
 #                                 last update: 3/28/07 {10:50:17 AM}
 # Author: Jonathan Guyer
 # E-mail: <guyer@nist.gov>
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
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 2007-03-28 JEG 1.0 original
 # 
 # #############################################################################
 ##

from fipy.terms.diffusionTerm import DiffusionTerm

class CollectedDiffusionTerm(DiffusionTerm):
    def __init__(self):
        self.orders = [None, None]

    def __iadd__(self, other):
        if not isinstance(other, DiffusionTerm):
            return DiffusionTerm.__iadd__(self, other)
        elif isinstance(other, CollectedDiffusionTerm):
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
        
    def _buildMatrix(self, var, boundaryConditions, dt, master=None):
        from fipy.tools import numerix
        from fipy.tools.sparseMatrix import _SparseMatrix
        
        N = len(var)
        RHSvector = numerix.zeros((N,),'d')
        matrix = _SparseMatrix(size=N)

        for term in self.orders:
            if term is not None:
                termMatrix, termRHSvector = term._buildMatrix(var, 
                                                              boundaryConditions, 
                                                              dt=dt, master=master)

                matrix += termMatrix
                RHSvector += termRHSvector
        
        return (matrix, RHSvector)

    def copy(self):
        dup = CollectedDiffusionTerm()
        dup.orders = self.orders
        return dup
        
    def _getDefaultSolver(self, solver):
        for term in self.orders:
            if term is not None:
                solver = term._getDefaultSolver(solver)
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
           DiffusionTerm(coeff=-1.0) + DiffusionTerm(coeff=(2.0, -3.0))

        """
        dup = CollectedDiffusionTerm()
        
        dup.orders = []
        for term in self.orders:
            if term is None:
                dup.orders.append(term)
            else:
                dup.orders.append(-term)
                
        return dup


