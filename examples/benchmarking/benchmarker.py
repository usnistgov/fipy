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
from tempfile import mkstemp
import popen2
import resource
        
from fipy.tests import doctestPlus
import examples.phase.anisotropy

N = 1000
steps = 5

script = doctestPlus._getScript("examples.phase.anisotropy")

script = script.replace("__main__", 
                        "__DONT_RUN_THIS__")
                        
script = script.replace("""nx = ny = 20""", 
                        """nx = ny = %d""" % N)

script = script.replace("""steps = 10""", 
                        """steps = 0""")

fd, path = mkstemp(".py")
os.write(fd, script)
os.close(fd)

p = popen2.Popen3('python "%s" --inline' % path)
p.wait()

ru0 = resource.getrusage(resource.RUSAGE_CHILDREN)

script = script.replace("""steps = 0""", 
                        """steps = %d""" % steps)

f = open(path, "w")
f.write(script)
f.close()

p = popen2.Popen4('python "%s"' % path)
p.wait()

for l in p.fromchild:
    print l.rstrip()

ru = resource.getrusage(resource.RUSAGE_CHILDREN)

print "-" * 79

print "          user time: %.9f s / step / cell" % ((ru.ru_utime - ru0.ru_utime) / steps / N**2)
print "        system time: %.9f s / step / cell" % ((ru.ru_stime - ru0.ru_utime) / steps / N**2)
print "max resident memory: %.2f B / cell" % (float(ru.ru_maxrss) / N**2)

os.remove(path)
