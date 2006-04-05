## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "meshIterator.py"
 #                                     created: 3/3/06 {9:00:00 PM}
 #                                 last update: 3/5/06 {6:34:06 PM}
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

class MeshIterator:
    def __init__(self, mesh, ids=(), checkIDs=False):
        self.mesh = mesh
        if type(ids) is type(1):
            ids = (ids,)
        ids = numerix.array(ids)
        if checkIDs and not self._canContain(ids):
            raise IndexError, 'Invalid IDs: %s' % str(other)
        self.ids = ids
            
    def __iter__(self):
        return iter(self.getIDs())
        
    def getMesh(self):
        return self.mesh
        
    def getIDs(self):
        return self.ids
        
    def __len__(self):
        return len(self.getIDs())
        
    def __array__(self, t = None):
        return numerix.array(self.getIDs(), t)
        
    def __repr__(self):
        return "%s(mesh=%s, ids=%s)" % (self.__class__.__name__,`self.getMesh()`, `self.getIDs()`)
        
    def __str__(self):
        return str(self.getIDs())
        
    def where(self, condition):
        return self.__class__(mesh=self.mesh, 
                              ids=numerix.compress(condition, self.getIDs()))
                              
    def __getitem__(self, index):
        return self.__class__(mesh=self.getMesh(), ids=self.getIDs()[index])
        
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
        return numerix.take(self.mesh.getFaceCenters(), self.getIDs())
        
    def getAreas(self):
        return numerix.take(self.mesh._getFaceAreas(), self.getIDs())

    def _getAdjacentCellIDs(self):
        return numerix.take(self.getMesh()._getAdjacentCellIDs()[0], self.getIDs())

##     def __init__(self, mesh, ids = ()):
##         self.mesh = mesh
##         self.ids = numerix.array(ids)
        
##     def __add__(self, other):
##         if isinstance(other, FaceIterator):
##             assert(self.mesh == other.getMesh())
##         return FaceIterator(mesh=self.mesh, 
##                             ids=self.ids + other)
                  
##     def __add__(self, other):
##         return FaceIterator(mesh=self.mesh, ids=list(self.ids) + list(other))

##     def __mul__(self, other):
##         return FaceIterator(mesh=self.mesh, ids=self.ids * other)
        
##     def __getitem__(self, index):
##         return self.ids[index]
        
##     def append(self, other):
##         ids = list(self.ids)
##         if isinstance(other,FaceIterator):
##             ids += list(other.ids)
##         else:
##             ids += list(other)
##             
##         return FaceIterator(mesh=self, ids=ids)
        

