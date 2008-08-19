#!/usr/bin/env python

import PyTrilinos
from PyTrilinos import Epetra
pid = str(Epetra.PyComm().MyPID())
from fipy import *
from fipy.meshes.numMesh.parallelGrid1D import ParallelGrid1D
import profiler


#fudge = profiler.calibrate_profiler(10000)

#profile = profiler.Profiler('profile'+str(pid), fudge=fudge)

pid = Epetra.PyComm().MyPID()

m = ParallelGrid1D(nx=10)

v = CellVariable(mesh=m,value=m.getCellCenters())

print repr(v)

e = TransientTerm() == DiffusionTerm()

e.solve(v,dt=.1)

print repr(v)
#profile.stop()

Epetra.Finalize()
