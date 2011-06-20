
import string
import time
import os
from copy_script import Copy_script

if os.path.exists('liquidVapor2D.py'):
    os.remove('liquidVapor2D.py')

DocProg = Copy_script(To='liquidVapor2D.py', From='examples/reactiveWetting/liquidVapor2D.py')
DocProg.finalize_options()
DocProg.run()
f = open('liquidVapor2D.py','r+w')
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
os.remove('liquidVapor2D.py')

g=open('liquidVapor2D.py', 'w')
g.write("".join(flist))
g.close()


