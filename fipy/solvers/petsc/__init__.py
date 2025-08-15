import logging

_log = logging.getLogger(__name__)

from .linearLUSolver import *
from .linearCGSolver import *
from .linearGMRESSolver import *
from .linearBicgSolver import *
from .linearCGSSolver import *
from .dummySolver import *
from . import petscConvergence

from .preconditioners import *

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = DummySolver
GeneralSolver = DefaultSolver
