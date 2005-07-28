
class Foo:
    def _bar(x,y,z): 'This is a sample method'

def baz(a,b):
    '''
    @param a: foo
    @type a: L{Foo}
    '''

from epydoc.html import HTMLFormatter
from epydoc.uid import *

import sys
print '~'*70
print '~'*70
print HTMLFormatter.format(sys)
print '~'*70
print '~'*70
#print HTMLFormatter.format(baz)
#print '~'*70
#print '~'*70
