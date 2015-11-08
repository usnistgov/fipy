#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "memoryLeak.py"
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
 
"""

This python script is ripped from
http://www.nightmare.com/medusa/memory-leaks.html

It outputs the top 100 number of outstanding references for each
object.

"""
from __future__ import print_function

__all__ = []

import sys
import types

def _get_refcounts(theClass = None):
    d = {}
    sys.modules
    # collect all classes
    for m in list(sys.modules.values()):
        for sym in dir(m):
            o = getattr (m, sym)
            if type(o) is type:
                if theClass is not None and o is not theClass:
                    continue
                d[o] = sys.getrefcount (o)
    # sort by refcount
    pairs = [(x[1],x[0]) for x in list(d.items())]
    pairs.sort()
    pairs.reverse()
    return pairs

def _print_top_N(n = 100, theClass = None):
    for n, c in _get_refcounts(theClass)[:n]:
        print('%10d %s' % (n, c.__name__))

if __name__ == '__main__':
    print_top_N()
