from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

from fipy.solvers.pyAMG.linearGMRESSolver import *
from fipy.solvers.pyAMG.linearCGSSolver import *
from fipy.solvers.pyAMG.linearPCGSolver import *
from fipy.solvers.pyAMG.linearLUSolver import *
from fipy.solvers.pyAMG.linearGeneralSolver import *

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = LinearGMRESSolver
GeneralSolver = LinearGeneralSolver
