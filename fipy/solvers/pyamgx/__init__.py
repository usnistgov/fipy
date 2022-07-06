from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

import atexit

import pyamgx

from .pyAMGXSolver import *
from .linearCGSolver import *
from .linearGMRESSolver import *
from .linearFGMRESSolver import *
from .linearBiCGStabSolver import *
from .linearLUSolver import *
from .aggregationAMGSolver import *
from .classicalAMGSolver import *

pyamgx.initialize()
atexit.register(pyamgx.finalize)

DefaultSolver = LinearCGSolver
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

__all__.extend(linearCGSolver.__all__)
__all__.extend(linearFGMRESSolver.__all__)
__all__.extend(linearBiCGStabSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(aggregationAMGSolver.__all__)
