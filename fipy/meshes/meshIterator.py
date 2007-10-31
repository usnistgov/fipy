## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "meshIterator.py"
 #                                     created: 3/3/06 {9:00:00 PM}
 #                                 last update: 10/27/07 {10:19:55 AM}
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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
 # 2006-03-04 JEG 1.0 original
 # 
 # ########################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix

class MeshIterator(list):
    def __init__(self, mesh, ids=(), checkIDs=False):
        if type(ids) is type(1) or numerix.shape(ids) is ():
            ids = (ids,)
        list.__init__(self, ids)
        self.mesh = mesh
        if checkIDs and not self._canContain(ids):
            raise IndexError, 'Invalid IDs: %s' % str(other)
        
    def getMesh(self):
        return self.mesh
        
    def getIDs(self):
        return self[:]
        
    def __repr__(self):
        return "%s(mesh=%s, ids=%s)" % (self.__class__.__name__,`self.getMesh()`, `self.getIDs()`)
        
    def __str__(self):
        return str(self.getIDs())
        
    def where(self, condition):
        return self.__class__(mesh=self.mesh, 
                              ids=numerix.compress(condition, self.getIDs()))
              
                              ##     def __getitem__(self, index):
                              ##         return self.__class__(mesh=self.getMesh(), ids=self.getIDs()[index])

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            other = self.__class__(mesh=self.getMesh(), ids=other)
            
        assert self.getMesh() == other.getMesh()

        return self.__class__(mesh=self.getMesh(), 
                              ids=list(self.getIDs()) + list(other.getIDs()))

class FaceIterator(MeshIterator):
    def _canContain(self, other):
        import sets
        return (sets.Set(other).issubset(sets.Set(self.getMesh().getFaces())))
        
    def getCenters(self):
        return numerix.take(self.mesh.getFaceCenters(), self.getIDs(), axis=1)
        
    def getAreas(self):
        return numerix.take(self.mesh._getFaceAreas(), self.getIDs())

    def _getAdjacentCellIDs(self):
        return numerix.take(self.getMesh()._getAdjacentCellIDs()[0], self.getIDs())
