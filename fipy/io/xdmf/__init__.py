from document import Open
from time import TimeSeries

def _xml_to_array(xml):
    from fipy.tools import numerix
    from StringIO import StringIO
    return numerix.loadtxt(StringIO(xml)).ravel()

