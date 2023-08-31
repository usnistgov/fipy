""":ref:`Discretizations <section:discretization>` of partial differential equation expressions
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

class ExplicitVariableError(Exception):
    def __init__(self, s='Terms with explicit Variables cannot mix with Terms with implicit Variables.'):
        Exception.__init__(self, s)

class TermMultiplyError(Exception):
    def __init__(self, s='Must multiply terms by int or float.'):
        Exception.__init__(self, s)

class AbstractBaseClassError(NotImplementedError):
    def __init__(self, s="can't instantiate abstract base class"):
        NotImplementedError.__init__(self, s)

class VectorCoeffError(TypeError):
    def __init__(self, s="The coefficient must be a vector value."):
        TypeError.__init__(self, s)

class SolutionVariableNumberError(Exception):
    def __init__(self, s='Different number of solution variables and equations.'):
        Exception.__init__(self, s)

class SolutionVariableRequiredError(Exception):
    def __init__(self, s='The solution variable needs to be specified.'):
        Exception.__init__(self, s)

class IncorrectSolutionVariable(Exception):
    def __init__(self, s='The solution variable is incorrect.'):
        Exception.__init__(self, s)

class TransientTermError(Exception):
    def __init__(self, s='The equation requires a TransientTerm with explicit convection.'):
        Exception.__init__(self, s)

from fipy.terms.transientTerm import *
from fipy.terms.diffusionTerm import *
from fipy.terms.explicitDiffusionTerm import *
from fipy.terms.implicitDiffusionTerm import *
from fipy.terms.implicitSourceTerm import *
from fipy.terms.residualTerm import *
from fipy.terms.centralDiffConvectionTerm import *
from fipy.terms.explicitUpwindConvectionTerm import *
from fipy.terms.exponentialConvectionTerm import *
from fipy.terms.hybridConvectionTerm import *
from fipy.terms.powerLawConvectionTerm import *
from fipy.terms.upwindConvectionTerm import *
from fipy.terms.vanLeerConvectionTerm import *
from fipy.terms.firstOrderAdvectionTerm import *
from fipy.terms.advectionTerm import *
ConvectionTerm = PowerLawConvectionTerm
