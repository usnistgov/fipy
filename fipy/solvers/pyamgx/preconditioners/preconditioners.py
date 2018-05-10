import copy

__all__ = ["AMGPreconditioner",
           "CGPreconditioner",
           "BiCGStabPreconditioner",
           "FGMRESPreconditioner",
           "BlockJacobiPreconditioner",
           "MultiColorDILUPreconditioner"]

class Preconditioner:
    def __init__(self, preconditioner_type):
        self.config_dict = {
            "solver": preconditioner_type,
            "max_iters": 1
        }
    def __call__(self, **kwargs):
        """
        :Parameters:
            - kwargs: Keyword arguments specifying AMGX solver options.
        """
        self.config_dict.update(kwargs)
        return copy.copy(self.config_dict)

AMGPreconditioner = Preconditioner("AMG")
CGPreconditioner = Preconditioner("PCG")
BiCGStabPreconditioner = Preconditioner("PCIBCGSTAB")
FGMRESPreconditioner = Preconditioner("FGMRES")
BlockJacobiPreconditioner = Preconditioner("BLOCK_JACOBI")
MultiColorDILUPreconditioner = Preconditioner("MULTICOLOR_DILU")
