from __future__ import print_function
from __future__ import unicode_literals
import os
import sys
import re
from subprocess import Popen, PIPE

from fipy import numerix

from fipy.tools.parser import parse

from examples.benchmarking.utils import monitor

def main():
    steps = parse('--numberOfSteps', action='store',
                  type='int', default=20)

    benchmarker = os.path.join(os.path.dirname(__file__),
                               "benchmarker.py")

    args = sys.argv[1:]

    print("size\tcpu / (s / step / cell)\trsz / (B / cell)\tvsz / (B / cell)")

    pyth = sys.executable or "python"

    for size in numerix.arange(2, 6.5, 0.5):
        p = Popen([pyth, benchmarker] + args
                  + ["--numberOfElements=%d" % int(10**size),
                     "--numberOfSteps=0"],
                  stdout=PIPE,
                  stderr=PIPE)

        cpu0, rsz0, vsz0 = monitor(p)

        p = Popen([pyth, benchmarker,
                   "--numberOfElements=%d" % int(10**size),
                   "--numberOfSteps=%d" % steps,
                   "--cpuBaseLine=%f" % cpu0] + args,
                  stdout=PIPE,
                  stderr=PIPE)

        cpu, rsz, vsz = monitor(p)

        print("%d\t%g\t%g\t%g" % (10**size, cpu, rsz, vsz))

if __name__ == "__main__":
    main()
