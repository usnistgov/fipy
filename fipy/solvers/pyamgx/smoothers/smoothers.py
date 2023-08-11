from __future__ import unicode_literals
from builtins import object
import copy

__all__ = ["BlockJacobiSmoother",
           "MultiColorDILUSmoother",
           "MultiColorGSSmoother",
           "MultiColorILUSmoother",
           "Smoother"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Smoother(object):
    """Interface to pyamgx_ `smoother configuration`_.

    .. _pyamgx: https://pyamgx.readthedocs.io
    .. _smoother configuration: https://pyamgx.readthedocs.io/en/latest/basic.html#config-objects
    """

    def __init__(self, smoother_type):
        self.config_dict = {
            "solver": smoother_type,
            "max_iters": 1
        }
    def __call__(self, **kwargs):
        self.config_dict.update(kwargs)
        return copy.copy(self.config_dict)

BlockJacobiSmoother = Smoother("BLOCK_JACOBI")
MultiColorDILUSmoother = Smoother("MULTICOLOR_DILU")
MultiColorGSSmoother = Smoother("MULTICOLOR_GS")
MultiColorILUSmoother = Smoother("MULTICOLOR_ILU")
