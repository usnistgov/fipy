#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
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

import popen2
import re

r, w = popen2.popen2("svn log --quiet --revision 3000:HEAD")

scanf_e = "[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?"

reUser = re.compile("user time: (%s) s / step / cell" % scanf_e)
reSystem = re.compile("system time: (%s) s / step / cell" % scanf_e)
reMemory = re.compile("max resident memory: (%s) B / cell" % scanf_e)

revs = re.findall(r"^r([0-9]*) \|", "".join(r), re.MULTILINE)

for rev in revs:
    r, w = popen2.popen2("svn update --revision %s" % rev)
    
    r, w = popen2.popen2("python benchmarker.py")
    
    r = "".join(r)
    
    user = reUser.search(r, re.MULTILINE)
    sys = reSystem.search(r, re.MULTILINE)
    mem = reMemory.search(r, re.MULTILINE)
    
    print rev, user.group(1), sys.group(1), mem.group(1)