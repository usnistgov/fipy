#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "sparseMatrix.py"
 #                                    created: 11/10/03 {3:15:38 PM} 
 #                                last update: 4/30/04 {3:40:42 PM} 
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


class SparseMatrix:
    
    """
    SparseMatrix class wrapper for pysparse.
    Allowing basic python operations __add__, __sub__ etx
    Facilitate matrix populating in an easy way
    An example subcalls will be the identity matrix or the empty matrix.
    """

    def __init__(self, size = None, vector = None, id1 = None, id2 = None, bandWidth = None, spmatrix = None):
        if spmatrix != None:
            self.matrix = spmatrix
        else
            self.matrix = spmatrix.ll_mat(size[0], size[1], size[0] * bandWidth)
            if vector != None:
                self.put(vector, id1, id2)
                
    def getMatrix():
        return self.matrix
    
    def copy(self):
	return self.matrix.copy()
	
    def __getitem__(self, index): 
	return self.matrix[index]
	
    def __str__(self):
	return self.matrix.__str__()
	    
    def __repr__(self):
	return self.matrix.__repr__()
	
    def __setitem__(self, index, value):
	self.matrix[index] = value
	
    def __add__(self, other):
        return self.matrix.copy().shift(1., other.getMatrix())

    def __sub__(self, other):
	return self.matrix.copy().shift(-1., other.getMatrix())

    def __rsub__(self, other):
        return other.matrix.copy().shift(-1., self.matrix)
	    
    def __mul__(self, other):
        if type(other) == type(self):
            return spmatrix.matrixMultiply(self.matrix, other.getMatrix())
        elif type(1) == type(other):
            
    def __neg__(self):
	return -self.matrix
	
    def __pos__(self):
	return +self.matrix
	
    def __eq__(self,other):
	return self.matrix.__eq__(other.getMatrix())

    def shape(self):
        return self.matrix.shape
	
    def transpose(self):
        return self.matrix.shape()

    def put(self, vector, id1, id2):
        if id2 == None:
            id2 = id1
        for i in range(len(vector)):
            self.matrix[id1[i], id2[i]] = vector[i]

    def take(self, id1, id2 = None):
        if id2 == None:
            id2 = id1
        vector = Numeric.zeros(len(id1), 'd')
        for i in range(len(id1)):
            vector[i] = self.matrix[id1[i], id2[i]]            
        return vector

    def putAdd(self, vector, id1, id2):
        tmp = SparseMatrix(matrix = 0. * self.copy())
        tmp.put(vector, id1, id2)
        return self + tmp
            
    def getNumeric(self):
        numMatrix = Numeric.zeros(self.getShape(), 'd')
        for i in self.shape[0]:
            for j in self.shape[1]:
                numMatrix[i, j] = self.matrix[i, j]

        return numMatrix

    def getDiagonal(self):
        numMatrix = Numeric.zeros(self.getShape()[0], 'd')
        for i in range(self.getShape[0])
