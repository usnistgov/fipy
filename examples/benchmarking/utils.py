from __future__ import unicode_literals
import re

import numpy

scanf_e = "[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?"

reCPU = re.compile("cpu time: (%s) s / step / cell" % scanf_e)
reRSZ = re.compile("max resident memory: (%s) B / cell" % scanf_e)
reVSZ = re.compile("max virtual memory: (%s) B / cell" % scanf_e)

def monitor(p):
    r = "".join(p.communicate()[0])

    cpu = reCPU.search(r, re.MULTILINE)
    rsz = reRSZ.search(r, re.MULTILINE)
    vsz = reVSZ.search(r, re.MULTILINE)

    def numOrNaN(m, g=1):
        if m is None:
            return numpy.NaN
        else:
            return float(m.group(g))

    return (numOrNaN(cpu),
            numOrNaN(rsz),
            numOrNaN(vsz))
