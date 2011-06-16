import string
import time
import os
from copy_script import Copy_script

if os.path.exists('mesh1D.py'):
    os.remove('mesh1D.py')

##call some function in copy_script class
DocProg = Copy_script(To='mesh1D.py', From='../../examples/diffusion/mesh1D.py')
DocProg.finalize_options()
DocProg.run()

f = open('mesh1D.py','r+w')
flist = f.readlines()
print flist

#for i in range(len(flist)):
for index, line in enumerate(flist):
    whitespaces = len(line) - len(line.lstrip())
    if 'from fipy import *' in line:
        flist.insert(index + 1, 'import time \ntimes = [] \ntimes += time.time()\n')
    if 'mesh =' in line and not '#' in line:
        flist.insert(index + 1, 'times += time.time()\n')
    if ('.sweep' in line or '.solve' in line) and not '#' in line:
        if whitespaces != 0:
            if '(' in line and not ')' in line:
                flist.insert(index + 2, whitespaces * ' ' + 'times += time.time()\n')
            else:
                flist.insert(index + 1, whitespaces * ' ' + 'times += time.time()\n')
        elif '(' in line and not ')' in line:
            flist.insert(index + 2, 'times += time.time()\n')
        else:
            flist.insert(index + 1, 'times += time.time()\n')
    if 'while' in line and not '#' in line:
        if whitespaces != 0:
            flist.insert(index + 1, whitespaces * ' ' + 'times += time.time()\n')
        else:
            flist.insert(index + 1, '    times += time.time()\n')
    if 'viewer.plot()' in line and not '#' in line:
       line_contents = line
       del line
       flist.insert(index, "##" + line_contents)
       

flist.append('\ntimes += time.time()\n')


f.close()
os.remove('mesh1D.py')

g=open('mesh1D.py', 'w')
g.write("".join(flist))
g.close()

##run benchmarks for fipy/tools/mesh1D.py, uploading results to Codespeed.


