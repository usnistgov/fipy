"""
#-----*-Pyth-*-
####################################################################
# PFM - Python-based phase field solver
#
# FILE: "variable.py"
#                                   created: 11/10/03 {3:13:28 PM} 
#                               last update: 11/13/03 {12:08:25 PM} 
# Author: Jonathan Guyer
# E-mail: guyer@nist.gov
#   mail: NIST
#    www: http://ctcms.nist.gov/
# 
#========================================================================
#This software was developed at the National Institute of Standards
#and Technology by employees of the Federal Government in the course
#of their official duties.  Pursuant to title 17 Section 105 of the
#United States Code this software is not subject to copyright
#protection and is in the public domain.  PFM is an experimental
#system.  NIST assumes no responsibility whatsoever for its use by
#other parties, and makes no guarantees, expressed or implied, about
#its quality, reliability, or any other characteristic.  We would
#appreciate acknowledgement if the software is used.
#
#This software can be redistributed and/or modified freely
#provided that any derivative works bear some notice that they are
#derived from it, and any modified versions bear some notice that
#they have been modified.
#========================================================================
# See the file "license.terms" for information on usage and  redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# 
####################################################################
#----
"""

import Numeric
import equation

class Variable:
    variables = ()
    
    def __init__(self,name,mesh,equation = equation.Equation(self,(,))):
	self.name = name
	self.mesh = mesh
	self.value = Numeric.zeroes([len(mesh.cells())],'d')
	self.equation = equation
	self.variables += self

