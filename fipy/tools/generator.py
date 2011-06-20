#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "serialCommWrapper.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Andrew Acquaviva <andrewa@nist.gov>
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

import string
import time
import os
from copy_script import Copy_script


##fromScripts = ['examples/cahnHilliard/mesh2D.py', 'examples/phase/anisotropy.py',\
##                   'examples/reactiveWetting/liquidVapor2D.py']
##toScripts= ['mesh2D.py', 'anisotropy.py', 'liquidVapor2D.py']

def run(extra_examples):
    fromScripts = ['examples/cahnHilliard/mesh2D.py', 'examples/phase/anisotropy.py',\
                   'examples/reactiveWetting/liquidVapor2D.py']
    toScripts= ['mesh2D.py', 'anisotropy.py', 'liquidVapor2D.py']
    fromScripts += extra_examples
    for i in extra_examples:
        deconstruct = i.split('/')
        toScripts += deconstruct[1+len(decsontruct)]
    for i in range(len(fromScripts)):
        if os.path.exists(toScripts[i]):
            os.remove(toScripts[i])
        print 'Running on file', fromScripts[i]
        DocProg = Copy_script(To=toScripts[i], From=fromScripts[i])
        DocProg.finalize_options()
        DocProg.run()
        f = open(toScripts[i],'r+w')
        flist = f.readlines()
    
    
        for index, line in enumerate(flist):
            whitespaces = len(line) - len(line.lstrip())
            if 'from fipy import *' in line:
                flist.insert(index + 1, 'import time \ntimes = [] \ntimes.append(time.time())\n')
            elif 'mesh =' in line and not '#' in line:
                flist.insert(index + 1, whitespaces * ' ' + 'times.append(time.time())\n')
            elif ('.sweep' in line or '.solve' in line) and not '#' in line:
                if whitespaces != 0:
                    if '(' in line and not ')' in line:
                        flist.insert(index + 2, whitespaces * ' ' + 'times.append(time.time())\n')
                    else:
                        flist.insert(index + 1, whitespaces * ' ' + 'times.append(time.time())\n')
                elif '(' in line and not ')' in line:
                    flist.insert(index + 2, 'times.append(time.time())\n')
                else:
                    flist.insert(index + 1, 'times.append(time.time())\n')
            elif 'while' in line and not '#' in line:
                if whitespaces != 0:
                    flist.insert(index + 1, (whitespaces + 4) * ' ' + 'times.append(time.time())\n')
                else:
                    flist.insert(index + 1, 4 * ' ' + 'times.append(time.time())\n')
            elif "__name__ == '__main__'" in line and not '#' in line:
                line = line.lstrip()
                split_line = line.lstrip('if ').split(':')
                commentedline = whitespaces * ' ' + 'if False and ' + split_line[0] + ":\n"
                flist.insert(index, commentedline)
                del flist[index+1]
            elif '__name__ == "__main__"' in line and not '#' in line:
                line = line.lstrip()
                split_line = line.lstrip('if ').split(':')
                commentedline = whitespaces * ' ' + 'if False and ' + split_line[0] + ":\n"
                flist.insert(index, commentedline)
                del flist[index+1]
        flist.append('\ntimes.append(time.time())\n')
        flist.append('\nruntime = times[len(times)-1]-times[0]')
        flist.append("\nprint 'runtime:', runtime") 

        f.close()
        os.remove(toScripts[i])
            
        g=open(toScripts[i], 'w')
        g.write("".join(flist))
        g.close()

if __name__ == "__main__":
    run([])
 
