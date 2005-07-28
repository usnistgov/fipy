#
# plaintext.py: plaintext docstring parsing
# Edward Loper
#
# Created [04/10/01 12:00 AM]
# $Id$
#

"""
Parser for plaintext docstrings.  Plaintext docstrings are rendered as
verbatim output, preserving all whitespace.
"""
from epydoc.markup import *

def parse_docstring(docstring, errors, **options):
    """
    @return: A pair C{(M{d}, M{e})}, where C{M{d}} is a
        C{ParsedDocstring} that encodes the contents of the given
        plaintext docstring; and C{M{e}} is a list of errors that were
        generated while parsing the docstring.
    @rtype: C{L{ParsedPlaintextDocstring}, C{list} of L{ParseError}}
    """
    return ParsedPlaintextDocstring(docstring, **options)

class ParsedPlaintextDocstring(ParsedDocstring):
    def __init__(self, text, **options):
        self._verbatim = options.get('verbatim', 1)
        self._text = text

    def to_html(self, docstring_linker, **options):
        if options.get('verbatim', self._verbatim) == 0:
            return plaintext_to_html(self.to_plaintext(docstring_linker))
        else:
            return ParsedDocstring.to_html(self, docstring_linker, **options)

    def to_latex(self, docstring_linker, **options):
        if options.get('verbatim', self._verbatim) == 0:
            return plaintext_to_latex(self.to_plaintext(docstring_linker))
        else:
            return ParsedDocstring.to_latex(self, docstring_linker, **options)

    def to_plaintext(self, docstring_linker, **options):
        return self._text
    
    def summary(self):
        m = re.match(r'(\s*[\w\W]*?\.)(\s|$)', self._text)
        if m:
            return ParsedPlaintextDocstring(m.group(1), verbatim=0)
        else:
            summary = self._text.split('\n', 1)[0]+'...'
            return ParsedPlaintextDocstring(summary, verbatim=0)
        
#     def concatenate(self, other):
#         if not isinstance(other, ParsedPlaintextDocstring):
#             raise ValueError, 'Could not concatenate docstrings'
#         text = self._text+other._text
#         options = self._options.copy()
#         options.update(other._options)
#         return ParsedPlaintextDocstring(text, options)
