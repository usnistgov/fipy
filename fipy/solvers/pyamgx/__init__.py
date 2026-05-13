
import logging

_log = logging.getLogger(__name__)

import atexit

import pyamgx

from .preconditioners import *

from .pyAMGXSolver import *
from .linearCGSolver import *
from .linearGMRESSolver import *
from .linearFGMRESSolver import *
from .linearBiCGStabSolver import *
from .linearLUSolver import *
from .aggregationAMGSolver import *
from .classicalAMGSolver import *

from . import pyamgxConvergence

pyamgx.initialize()
atexit.register(pyamgx.finalize)

DefaultSolver = LinearCGSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = DefaultSolver
GeneralSolver = LinearLUSolver
