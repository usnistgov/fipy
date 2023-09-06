from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from builtins import str
import os
import sys
import re
from tempfile import mkstemp
from subprocess import Popen, PIPE
from textwrap import dedent

from fipy.tools.parser import parse
from fipy.tools.numerix import sqrt
from fipy.tests import doctestPlus

import examples.phase.anisotropy

def main():
    numberOfElements = parse('--numberOfElements', action='store',
                             type='int', default=10000)
    N = int(sqrt(numberOfElements))

    steps = parse('--numberOfSteps', action='store',
                  type='int', default=20)

    start = parse('--startingStep', action='store',
                  type='int', default=0)

    cpu0 = parse('--cpuBaseLine', action='store',
                  type='float', default=0.)

    args = tuple(sys.argv[1:])

    dir = os.path.dirname(__file__)

    script = doctestPlus._getScript("examples.phase.anisotropy")

    script = script.replace("__main__",
                            "__DONT_RUN_THIS__")

    script = script.replace("nx = ny = 20",
                            "nx = ny = %d" % N)

    script = script.replace("steps = 10",
                            "steps = %d" % steps)

    fd, path = mkstemp(".py")
    os.write(fd, script)
    os.close(fd)

    cputime_RE = re.compile("(\d+:)?(\d+):(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?")

    def monitor(p):
        def cputimestring2secs(str):
            m = cputime_RE.match(str)

            secs = 60 * int(m.group(2)) + float(m.group(3))
            if m.group(1) is not None:
                secs += 3600 * int(m.group(1))

            return secs

        rsz = vsz = cpu = -1
        while p.poll() is None:
            try:
                r = Popen(("ps", "-p", str(p.pid), "-o" "rsz,vsz,cputime"),
                          stdout=PIPE,
                          env=os.environ).communicate()[0]

                ps = r.splitlines()[1].split()
                rsz = max(rsz, int(ps[0]) * 1024)
                vsz = max(vsz, int(ps[1]) * 1024)
                cpu = max(cpu, cputimestring2secs(ps[2]))
            except:
                break

        return rsz, vsz, cpu

    if start != 0:
        old = '''
              for i in range(steps):
              '''
        new = '''
              mesh_tmp, phase_tmp, dT_tmp = dump.read("%s/anisotropy-%d.dmp.gz")
              phase.setValue(phase_tmp.value)
              dT.setValue(dT_tmp.value)
              for i in range(steps):
              ''' % (dir, start)
        script = script.replace(dedent(old), dedent(new))

    old = '''\
          #    \cite{WarrenPolycrystal}.
          '''
    new = '''\
          #    \cite{WarrenPolycrystal}.

          dump.write((mesh, phase, dT), "%s/anisotropy-%%d.dmp.gz" %% (steps + %d))
          ''' % (dir, start)
    script = script.replace(dedent(old), dedent(new))

    f = open(path, "w")
    f.write(script)
    f.close()

    pyth = sys.executable or "python"

    p = Popen((pyth, path) + args, stdout=PIPE)

    rsz, vsz, cpu = monitor(p)

    for l in p.stdout:
        print(l.rstrip())

    print("-" * 79)

    if steps == 0:
        print("           cpu time: %.9f s / step / cell" % cpu)
        print("max resident memory: %.2f B / cell" % float(rsz))
        print(" max virtual memory: %.2f B / cell" % float(vsz))
    else:
        print("           cpu time: %.9f s / step / cell" % ((cpu - cpu0) / steps / N**2))
        print("max resident memory: %.2f B / cell" % (float(rsz) / N**2))
        print(" max virtual memory: %.2f B / cell" % (float(vsz) / N**2))

    os.remove(path)

if __name__ == "__main__":
    main()
