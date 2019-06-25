from future.utils import text_to_native_str
from future.utils import string_types

__all__ = [text_to_native_str("nativize_all")]

def nativize_all(t):
    def _nativize(s):
        if isinstance(s, string_types):
            s = text_to_native_str(s)
        return s
        
    return tuple([_nativize(s) for s in t])
    
