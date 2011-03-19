
import warnings

def getsetDeprecated(fn):
    """
    Used to decorate old get/set functions which are now deprecated in favor of
    their associated properties.
    """
    def new(*args, **kwargs):
        temp = "Get/set methods are deprecated; use the associated property instead."
        warnings.warn(temp, DeprecationWarning, stacklevel=3)

        return fn(*args, **kwargs)
    return new

