from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

from .preconditioners import *

from .linearCGSSolver import *
from .linearGMRESSolver import *
from .linearBicgstabSolver import *
from .linearLUSolver import *
from .linearPCGSolver import *

from . import scipyConvergence

DefaultSolver = LinearLUSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = LinearGMRESSolver
GeneralSolver = LinearLUSolver
