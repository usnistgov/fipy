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

from transientTerm import *
from diffusionTerm import *
from explicitDiffusionTerm import *
from implicitDiffusionTerm import *
from implicitSourceTerm import *
from residualTerm import *
from centralDiffConvectionTerm import *
from explicitUpwindConvectionTerm import *
from exponentialConvectionTerm import *
from hybridConvectionTerm import *
from powerLawConvectionTerm import *
from upwindConvectionTerm import *
from vanLeerConvectionTerm import *
ConvectionTerm = PowerLawConvectionTerm

__all__ = ["ExplicitVariableError",
           "TermMultiplyError",
           "AbstractBaseClassError",
           "VectorCoeffError",
           "SolutionVariableNumberError",
           "SolutionVariableRequiredError",
           "IncorrectSolutionVariable",
           "ConvectionTerm"]
           
__all__.extend(transientTerm.__all__)
__all__.extend(diffusionTerm.__all__)
__all__.extend(diffusionTermCorrection.__all__)
__all__.extend(diffusionTermNoCorrection.__all__)
__all__.extend(explicitDiffusionTerm.__all__)
__all__.extend(implicitDiffusionTerm.__all__)
__all__.extend(implicitSourceTerm.__all__)
__all__.extend(residualTerm.__all__)
__all__.extend(centralDiffConvectionTerm.__all__)
__all__.extend(explicitUpwindConvectionTerm.__all__)
__all__.extend(exponentialConvectionTerm.__all__)
__all__.extend(hybridConvectionTerm.__all__)
__all__.extend(powerLawConvectionTerm.__all__)
__all__.extend(upwindConvectionTerm.__all__)
__all__.extend(vanLeerConvectionTerm.__all__)
