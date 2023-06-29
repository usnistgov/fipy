from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

from fipy.solvers.scipy.linearCGSSolver import *
from fipy.solvers.scipy.linearGMRESSolver import *
from fipy.solvers.scipy.linearBicgstabSolver import *
from fipy.solvers.scipy.linearLUSolver import *
from fipy.solvers.scipy.linearPCGSolver import *

DefaultSolver = LinearLUSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = LinearGMRESSolver
GeneralSolver = LinearLUSolver
