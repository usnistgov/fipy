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
import shutil
from subprocess import Popen, PIPE
import sys
import tempfile

from fipy.tools.parser import parse

url = parse('--svnURL', action='store',
              type='string', default=None)

args = sys.argv[1:]

import pysvn

client = pysvn.Client()

# svn manipulations on the working copy in-place are dangerous

if url is None:
    info = client.info('.')
    url = info.url
    
dir = tempfile.mkdtemp()

env = os.environ.copy()
env['PYTHONPATH'] = dir

try:
    client.checkout(url, dir)

    scanf_e = "[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?"

    reCPU = re.compile("cpu time: (%s) s / step / cell" % scanf_e)
    reRSZ = re.compile("max resident memory: (%s) B / cell" % scanf_e)
    reVSZ = re.compile("max virtual memory: (%s) B / cell" % scanf_e)

    for entry in client.log(url):
        client.update(dir, revision=entry.revision)
        
        p = Popen(["python", 
                   os.path.join(os.path.dirname(__file__), 
                                "benchmarker.py")] + args, 
                  stdout=PIPE, 
                  stderr=PIPE, 
                  env=env)

        r = "".join(p.communicate()[0])
        
        cpu = reCPU.search(r, re.MULTILINE)
        rsz = reRSZ.search(r, re.MULTILINE)
        vsz = reVSZ.search(r, re.MULTILINE)
        
        print entry.revision.number, cpu.group(1), rsz.group(1), vsz.group(1)
except Exception, e:
    print e
    
shutil.rmtree(dir)
