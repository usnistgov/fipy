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
import re
from tempfile import mkstemp
import popen2
import resource
        
from fipy.tests import doctestPlus
import examples.phase.anisotropy

N = 100
steps = 20

script = doctestPlus._getScript("examples.phase.anisotropy")

script = script.replace("__main__", 
                        "__DONT_RUN_THIS__")
                        
script = script.replace("""nx = ny = 20""", 
                        """nx = ny = %d""" % N)

script = script.replace("""steps = 10""", 
                        """steps = 0""")

script0 = script.replace('''
#    \cite{WarrenPolycrystal}.''', '''
#    \cite{WarrenPolycrystal}.

dump.write((mesh, phase, dT), "anisotropy-0.dmp.gz")
''')

fd, path = mkstemp(".py")
os.write(fd, script0)
os.close(fd)

cputime_RE = re.compile("(\d+:)?(\d+):(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?")

def monitor(p):
    def cputimestring2secs(str):
        m = cputime_RE.match(str)
        
        secs = 60 * int(m.group(2)) + float(m.group(3))
        if m.group(1) is not None:
            secs += 3600 * int(m.group(1))
            
        return secs
        
    rsz = vsz = cputime = -1
    while p.poll() == -1:
        r, w = popen2.popen2("ps -p %d -o rsz,vsz,cputime" % p.pid)
        ps = r.readlines()[1].split()
        rsz = max(rsz, int(ps[0]) * 1024)
        vsz = max(vsz, int(ps[1]) * 1024)
        cputime = max(cputime, cputimestring2secs(ps[2]))
        
    return rsz, vsz, cputime
    
p = popen2.Popen3('python "%s" --inline' % path)

rsz0, vsz0, cputime0 = monitor(p)

script = script.replace("""steps = 0""", 
                        """steps = %d""" % steps)
                        
                        
datafile = file("data.txt", mode="w+", buffering=1)
datafile.write("step\cpu / (s / step / cell)\trsz / (B / cell)\tvsz / (B / cell)\n")

for block in range(500):
    script1 = script.replace('''
for i in range(steps):''', '''
mesh_tmp, phase_tmp, dT_tmp = dump.read("anisotropy-%d.dmp.gz")
phase.setValue(phase_tmp.getValue())
dT.setValue(dT_tmp.getValue())
for i in range(steps):''' % (block * steps))

    script1 = script1.replace('''
#    \cite{WarrenPolycrystal}.''', '''
#    \cite{WarrenPolycrystal}.
                    
dump.write((mesh, phase, dT), "anisotropy-%%d.dmp.gz" %% (steps + %d))
''' % (block * steps))

    f = open(path, "w")
    f.write(script1)
    f.close()

    p = popen2.Popen4('python "%s"' % path)
    
    rsz, vsz, cputime = monitor(p)

    for l in p.fromchild:
        print l.rstrip()

    print "-" * 79

    print "           cpu time: %.9f s / step / cell" % ((cputime - cputime0) / steps / N**2)
    print "max resident memory: %.2f B / cell" % (float(rsz) / N**2)
    print " max virtual memory: %.2f B / cell" % (float(vsz) / N**2)

    datafile.write("%d\t%g\t%g\t%g\n" % ((block + 1) * steps,
                                         (cputime - cputime0) / steps / N**2,
                                         float(rsz) / N**2,
                                         float(vsz) / N**2))
    
os.remove(path)

datafile.close()


