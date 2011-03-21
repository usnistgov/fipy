import re
import warnings

from numpy.lib.utils import _Deprecate

class _GetSetDeprecated(_Deprecate):
    def __init__(self, old_name=None, new_name=None, message=None):
        _Deprecate.__init__(self, old_name=old_name, new_name=None, message=message)
        self.new_property = new_name
        
    def __call__(self, func, *args, **kwargs):
        old_name = self.old_name
        new_name = self.new_name
        message = self.message

        if old_name is None:
            try:
                old_name = func.func_name
            except AttributeError:
                old_name = func.__name__
                
        self.old_name = old_name + "()"

        if new_name is None:
            RE = re.search("(_*)(get|set)(.)(.*)", old_name)
            if RE is not None:
                new_name = RE.group(1) + RE.group(3).lower() + RE.group(4)
                if message is None:
                    self.message = "Use the `%s` property instead." % new_name
        
        return _Deprecate.__call__(self, func=func, *args, **kwargs)

def getsetDeprecated(*args, **kwargs):
    if args:
        fn = args[0]
        args = args[1:]

        return _GetSetDeprecated(*args, **kwargs)(fn)
    else:
        return _GetSetDeprecated(*args, **kwargs)

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
