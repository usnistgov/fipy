import logging

_log = logging.getLogger(__name__)

from .linearLUSolver import *
from .linearPCGSolver import *
from .linearGMRESSolver import *
from .linearBicgSolver import *
from .linearCGSSolver import *
from .dummySolver import *

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = DummySolver
GeneralSolver = DefaultSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"]
           
__all__.extend(linearLUSolver.__all__)
__all__.extend(linearPCGSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearBicgSolver.__all__)
__all__.extend(linearCGSSolver.__all__)
