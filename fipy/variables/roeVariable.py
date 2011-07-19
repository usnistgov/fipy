#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "roeVariable.py"
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.tools import numerix
from fipy.variables.faceVariable import FaceVariable

class _RoeVariable(FaceVariable):
    def __init__(self, var, coeff):
        super(FaceVariable, self).__init__(mesh=var.mesh, elementshape=(2,) + var.shape[:-1], cached=True)
        self.var = var
##        self.var = self._requires(var)
        self.coeff = self._requires(coeff)
        
    def _calcValue(self):
        id1, id2 = self.mesh._adjacentCellIDs
        mesh = self.var.mesh
        
##        ## varDown.shape = (Nequ, Nfac)
##        varDown = numerix.take(self.var, id1, axis=-1)
##        varUp = numerix.take(self.var, id2, axis=-1)

        ## coeffDown.shape = (Nequ, Nequ, Nfac)
        coeffDown = (numerix.take(self.coeff, id1, axis=-1) * mesh._orientedAreaProjections[:, numerix.newaxis, numerix.newaxis]).sum(0)
        coeffUp = (numerix.take(self.coeff, id2, axis=-1) * mesh._orientedAreaProjections[:, numerix.newaxis, numerix.newaxis]).sum(0)
        
##         ## Anumerator.shape = (Nequ, Nface)
##         Anumerator = (coeffUp * varUp[:, numerix.newaxis] \
##                       - coeffDown * varDown[:, numerix.newaxis]).sum(1)

##         ## Adenominator.shape = (Nequ, Nfac)
##         Adenominator = varUp - varDown
##         Adenominator = numerix.where(Adenominator == 0,
##                                      1e-10,
##                                      Adenominator)

##         ## A.shape = (Nequ, Nequ, Nfac)
##         A = Anumerator[:, numerix.newaxis] / Adenominator[numerix.newaxis]

        ## Helper functions
        def mul(A, B):
            """
            Matrix multiply N MxM matrices, A.shape = B.shape = (M, M, N).
            """
            return numerix.sum(A.swapaxes(0,1)[:, :, numerix.newaxis] * B[:, numerix.newaxis], 0)

        def inv(A):
            """
            Inverts N MxM matrices, A.shape = (M, M, N).
            """
            return numerix.array(map(numerix.linalg.inv, A.transpose(2, 0, 1))).transpose(1, 2, 0)

        def eig(A):
            """
            Calculate the eigenvalues and eigenvectors of N MxM matrices, A.shape = (M, M, N).
            """
            tmp = zip(*map(numerix.linalg.eig, A.transpose(2, 0, 1)))
            return numerix.array(tmp[0]).swapaxes(0,1), numerix.array(tmp[1]).transpose(1,2,0)

        def sortedeig(A):
            """
            Caclulates the sorted eigenvalues and eigenvectors of N MxM matrices, A.shape = (M, M, N).
            """
            N = A.shape[-1]
            eigenvalues, R = eig(A)
            order = eigenvalues.argsort(0).swapaxes(0, 1)
            Nlist = [[i] for i in xrange(N)]
            return (eigenvalues[order, Nlist].swapaxes(0, 1),
                    R[:, order, Nlist].swapaxes(1, 2))

        ## A.shape = (Nequ, Nequ, Nfac)
        A = (coeffUp + coeffDown) / 2.
        
        eigenvalues, R = sortedeig(A)
        E = abs(eigenvalues) * numerix.identity(eigenvalues.shape[0])[..., numerix.newaxis]
        Abar = mul(mul(R, E), inv(R))

        ## value.shape = (2, Nequ, Nequ, Nfac)+
        value = numerix.zeros((2,) + A.shape, 'd')
        value[0] = (coeffDown + Abar) / 2
        value[1] = (coeffUp - Abar) / 2

        return value    

            
