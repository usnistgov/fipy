from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

from .linearCGSSolver import *
from .linearCGSolver import *
from .linearGMRESSolver import *
from .linearLUSolver import *
from .linearJORSolver import *

from .preconditioners import *

from . import pysparseConvergence

DefaultSolver = LinearCGSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = DefaultSolver
GeneralSolver =  LinearLUSolver
