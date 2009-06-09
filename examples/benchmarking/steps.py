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

import os
import sys
import re
from tempfile import mkstemp
import resource
from subprocess import Popen, PIPE
from textwrap import dedent
        

benchmarker = os.path.join(os.path.dirname(__file__), 
                           "benchmarker.py")

steps = 2

args = sys.argv[1:]

scanf_e = "[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?"

reCPU = re.compile("cpu time: (%s) s / step / cell" % scanf_e)
reRSZ = re.compile("max resident memory: (%s) B / cell" % scanf_e)
reVSZ = re.compile("max virtual memory: (%s) B / cell" % scanf_e)

def monitor(p):
    r = "".join(p.communicate()[0])

    cpu = reCPU.search(r, re.MULTILINE)
    rsz = reRSZ.search(r, re.MULTILINE)
    vsz = reVSZ.search(r, re.MULTILINE)

    return (float(cpu.group(1)),
            float(rsz.group(1)),
            float(vsz.group(1)))
    
p = Popen(["python", benchmarker, 
           "--numberOfSteps=0"] + args, 
          stdout=PIPE,
          stderr=PIPE)

cpu0, rsz0, vsz0 = monitor(p)

print "step\tcpu / (s / step / cell)\trsz / (B / cell)\tvsz / (B / cell)"

for block in range(10):
    p = Popen(["python", benchmarker, 
               "--numberOfSteps=%d" % steps,
               "--cpuBaseLine=%f" % cpu0] + args, 
              stdout=PIPE,
              stderr=PIPE)

    cpu, rsz, vsz = monitor(p)

    print "%d\t%g\t%g\t%g" % ((block + 1) * steps, cpu, rsz, vsz)


