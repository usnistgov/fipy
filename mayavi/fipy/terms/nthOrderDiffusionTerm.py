#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "nthOrderDiffusionTerm.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  FiPy
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm

class NthOrderDiffusionTerm(ImplicitDiffusionTerm):
    def __init__(self, coeff):
        """
        .. attention:: This class is deprecated. Use `ImplicitDiffusionTerm` instead.
        """
        import warnings
        warnings.warn("ImplicitDiffusionTerm should be used instead of NthOrderDiffusionTerm", DeprecationWarning, stacklevel=2)
        ImplicitDiffusionTerm.__init__(self, coeff = coeff)
        
class ExplicitNthOrderDiffusionTerm(ExplicitDiffusionTerm):
    def __init__(self, coeff):
        """
        .. attention:: This class is deprecated. Use `ExplicitDiffusionTerm` instead.
        """
        import warnings
        warnings.warn("ExplicitDiffusionTerm should be used instead of ExplicitNthOrderDiffusionTerm", DeprecationWarning, stacklevel=2)
        ExplicitDiffusionTerm.__init__(self, coeff = coeff)
    

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
