from __future__ import unicode_literals
from builtins import object
import copy

__all__ = ["JacobiPreconditioner",
           "DILUPreconditioner",
           "GaussSeidelPreconditioner",
           "ILUPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Preconditioner(object):
    def __init__(self, preconditioner_type):
        self.config_dict = {
            "solver": preconditioner_type,
            "max_iters": 1
        }
    def __call__(self, **kwargs):
        self.config_dict.update(kwargs)
        return copy.copy(self.config_dict)

JacobiPreconditioner = Preconditioner("BLOCK_JACOBI")
DILUPreconditioner = Preconditioner("MULTICOLOR_DILU")
GaussSeidelPreconditioner = Preconditioner("MULTICOLOR_GS")
ILUPreconditioner = Preconditioner("MULTICOLOR_ILU")
