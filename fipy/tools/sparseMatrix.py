#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "sparseMatrix.py"
 #                                    created: 11/10/03 {3:15:38 PM} 
 #                                last update: 5/17/04 {4:22:56 PM} 
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##


import Numeric

import spmatrix

class SparseMatrix:
    
    """
    SparseMatrix class wrapper for pysparse.
    SparseMatrix is always NxN
    Allowing basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way
    An example subcalls will be the identity matrix or the empty matrix.
    """

    def __init__(self, size = None, bandwidth = 0, matrix = None):
        if matrix != None:
            self.matrix = matrix
        else:
            self.matrix = spmatrix.ll_mat(size, size, size * bandwidth)
                
    def getMatrix(self):
        return self.matrix
    
    def copy(self):
	return SparseMatrix(matrix = self.matrix.copy())
	
    def __getitem__(self, index): 
	return self.matrix[index]
	
##     def __str__(self):
## 	return self.matrix.__str__()
	    
    def __repr__(self):
	return repr(self.getNumeric())
## 	return self.matrix.__repr__()
	
    def __setitem__(self, index, value):
	self.matrix[index] = value
	
    def __add__(self, other):
	L = self.matrix.copy()
	L.shift(1, other.getMatrix())
	return SparseMatrix(matrix = L)

    def __sub__(self, other):
	L = self.matrix.copy()
	L.shift(-1, other.getMatrix())
	return SparseMatrix(matrix = L)

    def __mul__(self, other):
        if type(other) == type(self):
            return SparseMatrix(matrix = spmatrix.matrixmultiply(self.matrix, other.getMatrix()))
        elif type(1) == type(other) or type(1.) == type(other):
	    N = self.matrix.shape[0]
	    L = spmatrix.ll_mat(N, N)
	    L.put(other * Numeric.ones(N))
	    return SparseMatrix(matrix = spmatrix.matrixmultiply(self.matrix, L))
	elif type(Numeric.ones(1)) == type(other):
	    y = other.copy()
	    self.matrix.matvec(other, y)
	    return y
	else:
 	    raise TypeError
            
    def __rmul__(self, other):
	if type(Numeric.ones(1)) == type(other):
	    y = other.copy()
	    self.matrix.matvec_transp(other, y)
	    return y
	else:
	    return self * other
	
    def __neg__(self):
	return self * -1
	
    def __pos__(self):
	return self
	
##     def __eq__(self,other):
## 	return self.matrix.__eq__(other.getMatrix())

    def getShape(self):
        return self.matrix.shape
	
    def transpose(self):
        return self.matrix.fuckEverythingUp()

    def put(self, vector, id1, id2):
	self.matrix.put(vector, id1, id2)

    def putDiagonal(self, vector):
	ids = Numeric.arange(len(vector))
	self.put(vector, ids, ids)

    def take(self, id1, id2):
	vector = Numeric.zeros(veclen, 'd')
	self.matrix.take(vector, id1, id2)
        return vector

    def takeDiagonal(self):
	ids = Numeric.arange(len(vector))
	return self.take(ids, ids)

    def addAt(self, vector, id1, id2):
	self.matrix.update_add_at(vector, id1, id2)

    def addAtDiagonal(self, vector):
	ids = Numeric.arange(len(vector))
	self.addAt(vector, ids, ids)
	
    def getNumeric(self):
	shape = self.getShape()
	indices = Numeric.indices(shape)
	numMatrix = self.take(indices[0].flat, indices[1].flat)
	
	return Numeric.reshape(numMatrix, shape)

class SparseIdentityMatrix(SparseMatrix):
    def __init__(self, size):
	SparseMatrix.__init__(self, size = size, bandwidth = 1)
	self.put(Numeric.ones(size))
    