#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "unaryTerm.py"
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

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.terms.term import Term

class _UnaryTerm(Term):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _getVars(self):
        return [self.var]

    def _getTransientVars(self):
        return []
                
    def _getUncoupledTerms(self):
        return [self]
    
    def __repr__(self):
        """
        The representation of a `Term` object is given by,
        
           >>> print __UnaryTerm(123.456)
           __UnaryTerm(coeff=123.456)

        """
        if self.var is None:
            varString = ''
        else:
            varString = ', var=%s' % repr(self.var)

        return "%s(coeff=%s%s)" % (self.__class__.__name__, repr(self.coeff), varString)

class __UnaryTerm(_UnaryTerm): 
    """
    Dummy subclass for tests
    """
    pass 

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
