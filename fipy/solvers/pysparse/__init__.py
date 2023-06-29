from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

from fipy.solvers.pysparse.linearCGSSolver import *
from fipy.solvers.pysparse.linearPCGSolver import *
from fipy.solvers.pysparse.linearGMRESSolver import *
from fipy.solvers.pysparse.linearLUSolver import *
from fipy.solvers.pysparse.linearJORSolver import *

from fipy.solvers.pysparse.preconditioners import *

DefaultSolver = LinearPCGSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = DefaultSolver
GeneralSolver =  LinearLUSolver
