from __future__ import unicode_literals

import logging

_log = logging.getLogger(__name__)

import atexit

import pyamgx

from fipy.solvers.pyamgx.pyAMGXSolver import *
from fipy.solvers.pyamgx.linearCGSolver import *
from fipy.solvers.pyamgx.linearGMRESSolver import *
from fipy.solvers.pyamgx.linearFGMRESSolver import *
from fipy.solvers.pyamgx.linearBiCGStabSolver import *
from fipy.solvers.pyamgx.linearLUSolver import *
from fipy.solvers.pyamgx.aggregationAMGSolver import *
from fipy.solvers.pyamgx.classicalAMGSolver import *

pyamgx.initialize()
atexit.register(pyamgx.finalize)

DefaultSolver = LinearCGSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = DefaultSolver
GeneralSolver = LinearLUSolver
