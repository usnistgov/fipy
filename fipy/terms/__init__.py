class ExplicitVariableError(Exception):
    def __init__(self):
        Exception.__init__(self, 'Terms with explicit Variables cannot mix with Terms with implicit Variables.')

class TermMultiplyError(Exception):
    def __init__(self):
        Exception.__init__(self, 'Must multiply terms by int or float."')

class AbstractBaseClassError(NotImplementedError):
    def __init__(self):
        NotImplementedError.__init__(self, "can't instantiate abstract base class")

class VectorCoeffError(TypeError):
    def __init__(self):
        TypeError.__init__(self, "The coefficient must be a vector value.")

class SolutionVariableNumberError(Exception):
    def __init__(self):
        Exception.__init__(self, 'Different number of solution variables and equations.')

class SolutionVariableRequiredError(Exception):
    def __init__(self):
        Exception.__init__(self, 'The solution variable needs to be specified.')

class IncorrectSolutionVariable(Exception):
    def __init__(self):
        Exception.__init__(self, 'The solution variable is incorrect.')

from transientTerm import TransientTerm

from diffusionTerm import DiffusionTerm
from diffusionTerm import DiffusionTermCorrection

from explicitDiffusionTerm import ExplicitDiffusionTerm
from implicitDiffusionTerm import ImplicitDiffusionTerm

from implicitSourceTerm import ImplicitSourceTerm

from residualTerm import ResidualTerm

from centralDiffConvectionTerm import CentralDifferenceConvectionTerm
from explicitUpwindConvectionTerm import ExplicitUpwindConvectionTerm
from exponentialConvectionTerm import ExponentialConvectionTerm
from hybridConvectionTerm import HybridConvectionTerm
from powerLawConvectionTerm import PowerLawConvectionTerm
from upwindConvectionTerm import UpwindConvectionTerm
from vanLeerConvectionTerm import VanLeerConvectionTerm
ConvectionTerm = PowerLawConvectionTerm

