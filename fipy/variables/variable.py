"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "variable.py"
                                   created: 11/10/03 {3:15:38 PM} 
                               last update: 11/20/03 {4:17:21 PM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
 Author: Daniel Wheeler
 E-mail: daniel.wheeler@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-10 JEG 1.0 original
###################################################################
"""

import Numeric

class Variable:
    
    def __init__(self, name, mesh, value=0.,viewer = 'None'):
	self.name = name
	self.mesh = mesh
	self.array = Numeric.zeros([len(mesh.getCells())],'d')
        self.viewer = viewer
        if viewer != 'None':
            self.viewer.setVar(self)
	    
	self.setValue(value)

    def plot(self):
        self.viewer.plot()
        
    def getMesh(self):
        return self.mesh

    def getArray(self):
        return self.array

    def getGridArray(self):
        return self.mesh.makeGridData(self.array)
	
    def setValue(self,value,cells = ()):
	if cells == ():
	    cells = self.mesh.getCells()
	    



    
    
