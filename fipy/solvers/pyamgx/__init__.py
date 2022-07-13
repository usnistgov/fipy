from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

import atexit

import pyamgx

from .preconditioners import *

from .pyAMGXSolver import *
from .linearPCGSolver import *
from .linearGMRESSolver import *
from .linearFGMRESSolver import *
from .linearBiCGStabSolver import *
from .linearLUSolver import *
from .aggregationAMGSolver import *
from .classicalAMGSolver import *

from . import pyamgxConvergence

pyamgx.initialize()
atexit.register(pyamgx.finalize)

DefaultSolver = LinearPCGSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = DefaultSolver
GeneralSolver = LinearLUSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"
          ]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

__all__.extend(linearPCGSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearFGMRESSolver.__all__)
__all__.extend(linearBiCGStabSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(aggregationAMGSolver.__all__)
__all__.extend(classicalAMGSolver.__all__)
__all__.extend(preconditioners.__all__)
