import copy

__all__ = ["BlockJacobiSmoother",
           "MultiColorDILUSmoother",
           "MultiColorGSSmoother",
           "MultiColorILUSmoother"]

class Smoother:
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
