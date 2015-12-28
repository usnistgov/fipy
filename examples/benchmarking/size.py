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
from subprocess import Popen, PIPE

from fipy import numerix

from fipy.tools.parser import parse

from examples.benchmarking.utils import monitor

steps = parse('--numberOfSteps', action='store',
              type='int', default=20)

benchmarker = os.path.join(os.path.dirname(__file__),
                           "benchmarker.py")

args = sys.argv[1:]

print "size\tcpu / (s / step / cell)\trsz / (B / cell)\tvsz / (B / cell)"

for size in numerix.arange(2,6.5,0.5):
    p = Popen(["python", benchmarker] + args
              + ["--numberOfElements=%d" % int(10**size),
                 "--numberOfSteps=0"],
              stdout=PIPE,
              stderr=PIPE)

    cpu0, rsz0, vsz0 = monitor(p)

    p = Popen(["python", benchmarker,
               "--numberOfElements=%d" % int(10**size),
               "--numberOfSteps=%d" % steps,
               "--cpuBaseLine=%f" % cpu0] + args,
              stdout=PIPE,
              stderr=PIPE)

    cpu, rsz, vsz = monitor(p)

    print "%d\t%g\t%g\t%g" % (10**size, cpu, rsz, vsz)
