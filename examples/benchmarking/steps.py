#!/usr/bin/env python


import os
import sys
import re
from subprocess import Popen, PIPE

from fipy.tools.parser import parse

from examples.benchmarking.utils import monitor

steps = parse('--numberOfSteps', action='store',
              type='int', default=20)

blocks = parse('--numberOfBlocks', action='store',
               type='int', default=20)

benchmarker = os.path.join(os.path.dirname(__file__),
                           "benchmarker.py")

args = sys.argv[1:]

p = Popen(["python", benchmarker] + args
          + ["--numberOfSteps=0"],
          stdout=PIPE,
          stderr=PIPE)

cpu0, rsz0, vsz0 = monitor(p)

print "step\tcpu / (s / step / cell)\trsz / (B / cell)\tvsz / (B / cell)"

for block in range(blocks):
    p = Popen(["python", benchmarker,
               "--startingStep=%d" % (block * steps),
               "--cpuBaseLine=%f" % cpu0] + args,
              stdout=PIPE,
              stderr=PIPE)

    cpu, rsz, vsz = monitor(p)

    print "%d\t%g\t%g\t%g" % (block * steps, cpu, rsz, vsz)
