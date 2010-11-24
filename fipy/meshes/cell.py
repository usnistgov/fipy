#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cell.py"
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

class Cell:
    def __init__(self, mesh, id):
        self.id = id
        self.mesh = mesh

    def getID(self):
        return self.id

    def getCenter(self):
        return self.mesh.getCellCenters()[...,self.id]

    def _getCellToCellDistances(self):
        return self.mesh._getCellToCellDistances()[...,self.id]

    def _getCellToCellIDs(self):
        return self.mesh._getCellToCellIDs()[...,self.id]

    def __cmp__(self, cell):
        return cmp(self.id, cell.getID())

    def getMesh(self):
        return self.mesh

    def getNormal(self, index):
        dis = self._getCellToCellDistances()[...,index]
        adjCellID = self._getCellToCellIDs()[...,index]
        vec = self.getCenter() - self.mesh.getCellCenters()[...,adjCellID]
        return vec / dis

    def __repr__(self):
        return "%s(mesh=%s, id=%s)" % (self.__class__.__name__,`self.getMesh()`, `self.getID()`)

                
                
