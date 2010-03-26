"""LaTeX.py

Character translation utilities for LaTeX-formatted text.

Usage:
 - unicode(string,'latex')
 - ustring.decode('latex')
are both available just by letting "import latex" find this file.
 - unicode(string,'latex+latin1')
 - ustring.decode('latex+latin1')
where latin1 can be replaced by any other known encoding, also
become available by calling latex.register().

We also make public a dictionary latex_equivalents,
mapping ord(unicode char) to LaTeX code.


Copyright (c) 2003,2008 David Eppstein

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
"""

import codecs
import re

def register():
    """Enable encodings of the form 'latex+x' where x describes another encoding.
    Unicode characters are translated to or from x when possible, otherwise
    expanded to latex.
    """
    codecs.register(_registry)

def getregentry():
    """Encodings module API."""
    return _registry('latex')

def _registry(encoding):
    if encoding == 'latex':
        encoding = None
    elif encoding.startswith('latex+'):
        encoding = encoding[6:]
    else:
        return None
        
    class Codec(codecs.Codec):
        def encode(self,input,errors='strict'):
            """Convert unicode string to latex."""
            output = []
            for c in input:
                if encoding:
                    try:
                        output.append(c.encode(encoding))
                        continue
                    except:
                        pass
                if ord(c) in latex_equivalents:
                    output.append(latex_equivalents[ord(c)])
                else:
                    output += ['{\\char', str(ord(c)), '}']
            return ''.join(output), len(input)
            
        def decode(self,input,errors='strict'):
            """Convert latex source string to unicode."""
            if encoding:
                input = unicode(input,encoding,errors)

            # Note: we may get buffer objects here.
            # It is not permussable to call join on buffer objects
            # but we can make them joinable by calling unicode.
            # This should always be safe since we are supposed
            # to be producing unicode output anyway.
            x = map(unicode,_unlatex(input))
            return u''.join(x), len(input)
    
    class StreamWriter(Codec,codecs.StreamWriter):
        pass
            
    class StreamReader(Codec,codecs.StreamReader):
        pass

    return (Codec().encode,Codec().decode,StreamReader,StreamWriter)

def _tokenize(tex):
    """Convert latex source into sequence of single-token substrings."""
    start = 0
    try:
        # skip quickly across boring stuff
        pos = _stoppers.finditer(tex).next().span()[0]
    except StopIteration:
        yield tex
        return

    while 1:
        if pos > start:
            yield tex[start:pos]
            if tex[start] == '\\' and not (tex[pos-1].isdigit() and tex[start+1].isalpha()):
                while pos < len(tex) and tex[pos].isspace(): # skip blanks after csname
                    pos += 1

        while pos < len(tex) and tex[pos] in _ignore:
            pos += 1    # flush control characters
        if pos >= len(tex):
            return
        start = pos
        if tex[pos:pos+2] in {'$$':None, '/~':None}:    # protect ~ in urls
            pos += 2
        elif tex[pos].isdigit():
            while pos < len(tex) and tex[pos].isdigit():
                pos += 1
        elif tex[pos] == '-':
            while pos < len(tex) and tex[pos] == '-':
                pos += 1
        elif tex[pos] != '\\' or pos == len(tex) - 1:
            pos += 1
        elif not tex[pos+1].isalpha():
            pos += 2
        else:
            pos += 1
            while pos < len(tex) and tex[pos].isalpha():
                pos += 1
            if tex[start:pos] == '\\char' or tex[start:pos] == '\\accent':
                while pos < len(tex) and tex[pos].isdigit():
                    pos += 1

class _unlatex:
    """Convert tokenized tex into sequence of unicode strings.  Helper for decode()."""

    def __iter__(self):
        """Turn self into an iterator.  It already is one, nothing to do."""
        return self

    def __init__(self,tex):
        """Create a new token converter from a string."""
        self.tex = tuple(_tokenize(tex))  # turn tokens into indexable list
        self.pos = 0                    # index of first unprocessed token 
        self.lastoutput = 'x'           # lastoutput must always be nonempty string
    
    def __getitem__(self,n):
        """Return token at offset n from current pos."""
        p = self.pos + n
        t = self.tex
        return p < len(t) and t[p] or None

    def next(self):
        """Find and return another piece of converted output."""
        if self.pos >= len(self.tex):
            raise StopIteration
        nextoutput = self.chunk()
        if self.lastoutput[0] == '\\' and self.lastoutput[-1].isalpha() and nextoutput[0].isalpha():
            nextoutput = ' ' + nextoutput   # add extra space to terminate csname
        self.lastoutput = nextoutput
        return nextoutput
    
    def chunk(self):
        """Grab another set of input tokens and convert them to an output string."""
        for delta,c in self.candidates(0):
            if c in _l2u:
                self.pos += delta
                return unichr(_l2u[c])
            elif len(c) == 2 and c[1] == 'i' and (c[0],'\\i') in _l2u:
                self.pos += delta       # correct failure to undot i
                return unichr(_l2u[(c[0],'\\i')])
            elif len(c) == 1 and c[0].startswith('\\char') and c[0][5:].isdigit():
                self.pos += delta
                return unichr(int(c[0][5:]))
    
        # nothing matches, just pass through token as-is
        self.pos += 1
        return self[-1]
    
    def candidates(self,offset):
        """Generate pairs delta,c where c is a token or tuple of tokens from tex
        (after deleting extraneous brackets starting at pos) and delta
        is the length of the tokens prior to bracket deletion.
        """
        t = self[offset]
        if t in _blacklist:
            return
        elif t == '{':
            for delta,c in self.candidates(offset+1):
                if self[offset+delta+1] == '}':
                    yield delta+2,c
        elif t == '\\mbox':
            for delta,c in self.candidates(offset+1):
                yield delta+1,c
        elif t == '$' and self[offset+2] == '$':
            yield 3, (t,self[offset+1],t)
        else:
            q = self[offset+1]
            if q == '{' and self[offset+3] == '}':
                yield 4, (t,self[offset+2])
            elif q:
                yield 2, (t,q)
            yield 1, t

latex_equivalents = {
    9: ' ',
    ord(u'\N{EN DASH}'): '{--}',
    ord(u'\N{EM DASH}'): '{---}',
    ord(u'\N{LEFT SINGLE QUOTATION MARK}'): '{`}',
    ord(u'\N{RIGHT SINGLE QUOTATION MARK}'): "{'}",
    ord(u'\N{LEFT DOUBLE QUOTATION MARK}'): '{``}',
    ord(u'\N{RIGHT DOUBLE QUOTATION MARK}'): "{''}",
    ord(u'\N{DAGGER}'): '{\\dag}',
    ord(u'\N{DOUBLE DAGGER}'): '{\\ddag}',
    ord(u'\N{BULLET}'): '{\\mbox{$\\bullet$}}',
    ord(u'\N{NUMBER SIGN}'): '{\\#}',
    ord(u'\N{AMPERSAND}'): '{\\&}',
    ord(u'\N{NO-BREAK SPACE}'): '{~}',
    ord(u'\N{INVERTED EXCLAMATION MARK}'): '{!`}',
    ord(u'\N{CENT SIGN}'): '{\\not{c}}',
    ord(u'\N{POUND SIGN}'): '{\\pounds}',
    ord(u'\N{SECTION SIGN}'): '{\\S}',
    ord(u'\N{DIAERESIS}'): '{\\"{}}',
    ord(u'\N{NOT SIGN}'): '{\\neg}',
    ord(u'\N{SOFT HYPHEN}'): '{\\-}',
    ord(u'\N{MACRON}'): '{\\={}}',
    ord(u'\N{DEGREE SIGN}'): '{\\mbox{$^\\circ$}}',
    ord(u'\N{PLUS-MINUS SIGN}'): '{\\mbox{$\\pm$}}',
    ord(u'\N{SUPERSCRIPT TWO}'): '{\\mbox{$^2$}}',
    ord(u'\N{SUPERSCRIPT THREE}'): '{\\mbox{$^3$}}',
    ord(u'\N{ACUTE ACCENT}'): "{\\'{}}",
    ord(u'\N{MICRO SIGN}'): '{\\mbox{$\\mu$}}',
    ord(u'\N{PILCROW SIGN}'): '{\\P}',
    ord(u'\N{MIDDLE DOT}'): '{\\mbox{$\\cdot$}}',
    ord(u'\N{CEDILLA}'): '{\\c{}}',
    ord(u'\N{SUPERSCRIPT ONE}'): '{\\mbox{$^1$}}',
    ord(u'\N{INVERTED QUESTION MARK}'): '{?`}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH GRAVE}'): '{\\`A}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH CIRCUMFLEX}'): '{\\^A}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH TILDE}'): '{\\~A}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH DIAERESIS}'): '{\\"A}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH RING ABOVE}'): '{\\AA}',
    ord(u'\N{LATIN CAPITAL LETTER AE}'): '{\\AE}',
    ord(u'\N{LATIN CAPITAL LETTER C WITH CEDILLA}'): '{\\c{C}}',
    ord(u'\N{LATIN CAPITAL LETTER E WITH GRAVE}'): '{\\`E}',
    ord(u'\N{LATIN CAPITAL LETTER E WITH ACUTE}'): "{\\'E}",
    ord(u'\N{LATIN CAPITAL LETTER E WITH CIRCUMFLEX}'): '{\\^E}',
    ord(u'\N{LATIN CAPITAL LETTER E WITH DIAERESIS}'): '{\\"E}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH GRAVE}'): '{\\`I}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH CIRCUMFLEX}'): '{\\^I}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH DIAERESIS}'): '{\\"I}',
    ord(u'\N{LATIN CAPITAL LETTER N WITH TILDE}'): '{\\~N}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH GRAVE}'): '{\\`O}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH ACUTE}'): "{\\'O}",
    ord(u'\N{LATIN CAPITAL LETTER O WITH CIRCUMFLEX}'): '{\\^O}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH TILDE}'): '{\\~O}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH DIAERESIS}'): '{\\"O}',
    ord(u'\N{MULTIPLICATION SIGN}'): '{\\mbox{$\\times$}}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH STROKE}'): '{\\O}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH GRAVE}'): '{\\`U}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH ACUTE}'): "{\\'U}",
    ord(u'\N{LATIN CAPITAL LETTER U WITH CIRCUMFLEX}'): '{\\^U}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH DIAERESIS}'): '{\\"U}',
    ord(u'\N{LATIN CAPITAL LETTER Y WITH ACUTE}'): "{\\'Y}",
    ord(u'\N{LATIN SMALL LETTER SHARP S}'): '{\\ss}',
    ord(u'\N{LATIN SMALL LETTER A WITH GRAVE}'): '{\\`a}',
    ord(u'\N{LATIN SMALL LETTER A WITH ACUTE}'): "{\\'a}",
    ord(u'\N{LATIN SMALL LETTER A WITH CIRCUMFLEX}'): '{\\^a}',
    ord(u'\N{LATIN SMALL LETTER A WITH TILDE}'): '{\\~a}',
    ord(u'\N{LATIN SMALL LETTER A WITH DIAERESIS}'): '{\\"a}',
    ord(u'\N{LATIN SMALL LETTER A WITH RING ABOVE}'): '{\\aa}',
    ord(u'\N{LATIN SMALL LETTER AE}'): '{\\ae}',
    ord(u'\N{LATIN SMALL LETTER C WITH CEDILLA}'): '{\\c{c}}',
    ord(u'\N{LATIN SMALL LETTER E WITH GRAVE}'): '{\\`e}',
    ord(u'\N{LATIN SMALL LETTER E WITH ACUTE}'): "{\\'e}",
    ord(u'\N{LATIN SMALL LETTER E WITH CIRCUMFLEX}'): '{\\^e}',
    ord(u'\N{LATIN SMALL LETTER E WITH DIAERESIS}'): '{\\"e}',
    ord(u'\N{LATIN SMALL LETTER I WITH GRAVE}'): '{\\`\\i}',
    ord(u'\N{LATIN SMALL LETTER I WITH ACUTE}'): "{\\'\\i}",
    ord(u'\N{LATIN SMALL LETTER I WITH CIRCUMFLEX}'): '{\\^\\i}',
    ord(u'\N{LATIN SMALL LETTER I WITH DIAERESIS}'): '{\\"\\i}',
    ord(u'\N{LATIN SMALL LETTER N WITH TILDE}'): '{\\~n}',
    ord(u'\N{LATIN SMALL LETTER O WITH GRAVE}'): '{\\`o}',
    ord(u'\N{LATIN SMALL LETTER O WITH ACUTE}'): "{\\'o}",
    ord(u'\N{LATIN SMALL LETTER O WITH CIRCUMFLEX}'): '{\\^o}',
    ord(u'\N{LATIN SMALL LETTER O WITH TILDE}'): '{\\~o}',
    ord(u'\N{LATIN SMALL LETTER O WITH DIAERESIS}'): '{\\"o}',
    ord(u'\N{DIVISION SIGN}'): '{\\mbox{$\\div$}}',
    ord(u'\N{LATIN SMALL LETTER O WITH STROKE}'): '{\\o}',
    ord(u'\N{LATIN SMALL LETTER U WITH GRAVE}'): '{\\`u}',
    ord(u'\N{LATIN SMALL LETTER U WITH ACUTE}'): "{\\'u}",
    ord(u'\N{LATIN SMALL LETTER U WITH CIRCUMFLEX}'): '{\\^u}',
    ord(u'\N{LATIN SMALL LETTER U WITH DIAERESIS}'): '{\\"u}',
    ord(u'\N{LATIN SMALL LETTER Y WITH ACUTE}'): "{\\'y}",
    ord(u'\N{LATIN SMALL LETTER Y WITH DIAERESIS}'): '{\\"y}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH MACRON}'): '{\\=A}',
    ord(u'\N{LATIN SMALL LETTER A WITH MACRON}'): '{\\=a}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH BREVE}'): '{\\u{A}}',
    ord(u'\N{LATIN SMALL LETTER A WITH BREVE}'): '{\\u{a}}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH OGONEK}'): '{\\c{A}}',
    ord(u'\N{LATIN SMALL LETTER A WITH OGONEK}'): '{\\c{a}}',
    ord(u'\N{LATIN CAPITAL LETTER C WITH ACUTE}'): "{\\'C}",
    ord(u'\N{LATIN SMALL LETTER C WITH ACUTE}'): "{\\'c}",
    ord(u'\N{LATIN CAPITAL LETTER C WITH CIRCUMFLEX}'): '{\\^C}',
    ord(u'\N{LATIN SMALL LETTER C WITH CIRCUMFLEX}'): '{\\^c}',
    ord(u'\N{LATIN CAPITAL LETTER C WITH DOT ABOVE}'): '{\\.C}',
    ord(u'\N{LATIN SMALL LETTER C WITH DOT ABOVE}'): '{\\.c}',
    ord(u'\N{LATIN CAPITAL LETTER C WITH CARON}'): '{\\v{C}}',
    ord(u'\N{LATIN SMALL LETTER C WITH CARON}'): '{\\v{c}}',
    ord(u'\N{LATIN CAPITAL LETTER D WITH CARON}'): '{\\v{D}}',
    ord(u'\N{LATIN SMALL LETTER D WITH CARON}'): '{\\v{d}}',
    ord(u'\N{LATIN CAPITAL LETTER E WITH MACRON}'): '{\\=E}',
    ord(u'\N{LATIN SMALL LETTER E WITH MACRON}'): '{\\=e}',
    ord(u'\N{LATIN CAPITAL LETTER E WITH BREVE}'): '{\\u{E}}',
    ord(u'\N{LATIN SMALL LETTER E WITH BREVE}'): '{\\u{e}}',
    ord(u'\N{LATIN CAPITAL LETTER E WITH DOT ABOVE}'): '{\\.E}',
    ord(u'\N{LATIN SMALL LETTER E WITH DOT ABOVE}'): '{\\.e}',
    ord(u'\N{LATIN CAPITAL LETTER E WITH OGONEK}'): '{\\c{E}}',
    ord(u'\N{LATIN SMALL LETTER E WITH OGONEK}'): '{\\c{e}}',
    ord(u'\N{LATIN CAPITAL LETTER E WITH CARON}'): '{\\v{E}}',
    ord(u'\N{LATIN SMALL LETTER E WITH CARON}'): '{\\v{e}}',
    ord(u'\N{LATIN CAPITAL LETTER G WITH CIRCUMFLEX}'): '{\\^G}',
    ord(u'\N{LATIN SMALL LETTER G WITH CIRCUMFLEX}'): '{\\^g}',
    ord(u'\N{LATIN CAPITAL LETTER G WITH BREVE}'): '{\\u{G}}',
    ord(u'\N{LATIN SMALL LETTER G WITH BREVE}'): '{\\u{g}}',
    ord(u'\N{LATIN CAPITAL LETTER G WITH DOT ABOVE}'): '{\\.G}',
    ord(u'\N{LATIN SMALL LETTER G WITH DOT ABOVE}'): '{\\.g}',
    ord(u'\N{LATIN CAPITAL LETTER G WITH CEDILLA}'): '{\\c{G}}',
    ord(u'\N{LATIN SMALL LETTER G WITH CEDILLA}'): '{\\c{g}}',
    ord(u'\N{LATIN CAPITAL LETTER H WITH CIRCUMFLEX}'): '{\\^H}',
    ord(u'\N{LATIN SMALL LETTER H WITH CIRCUMFLEX}'): '{\\^h}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH TILDE}'): '{\\~I}',
    ord(u'\N{LATIN SMALL LETTER I WITH TILDE}'): '{\\~\\i}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH MACRON}'): '{\\=I}',
    ord(u'\N{LATIN SMALL LETTER I WITH MACRON}'): '{\\=\\i}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH BREVE}'): '{\\u{I}}',
    ord(u'\N{LATIN SMALL LETTER I WITH BREVE}'): '{\\u\\i}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH OGONEK}'): '{\\c{I}}',
    ord(u'\N{LATIN SMALL LETTER I WITH OGONEK}'): '{\\c{i}}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH DOT ABOVE}'): '{\\.I}',
    ord(u'\N{LATIN SMALL LETTER DOTLESS I}'): '{\\i}',
    ord(u'\N{LATIN CAPITAL LIGATURE IJ}'): '{IJ}',
    ord(u'\N{LATIN SMALL LIGATURE IJ}'): '{ij}',
    ord(u'\N{LATIN CAPITAL LETTER J WITH CIRCUMFLEX}'): '{\\^J}',
    ord(u'\N{LATIN SMALL LETTER J WITH CIRCUMFLEX}'): '{\\^\\j}',
    ord(u'\N{LATIN CAPITAL LETTER K WITH CEDILLA}'): '{\\c{K}}',
    ord(u'\N{LATIN SMALL LETTER K WITH CEDILLA}'): '{\\c{k}}',
    ord(u'\N{LATIN CAPITAL LETTER L WITH ACUTE}'): "{\\'L}",
    ord(u'\N{LATIN SMALL LETTER L WITH ACUTE}'): "{\\'l}",
    ord(u'\N{LATIN CAPITAL LETTER L WITH CEDILLA}'): '{\\c{L}}',
    ord(u'\N{LATIN SMALL LETTER L WITH CEDILLA}'): '{\\c{l}}',
    ord(u'\N{LATIN CAPITAL LETTER L WITH CARON}'): '{\\v{L}}',
    ord(u'\N{LATIN SMALL LETTER L WITH CARON}'): '{\\v{l}}',
    ord(u'\N{LATIN CAPITAL LETTER L WITH STROKE}'): '{\\L}',
    ord(u'\N{LATIN SMALL LETTER L WITH STROKE}'): '{\\l}',
    ord(u'\N{LATIN CAPITAL LETTER N WITH ACUTE}'): "{\\'N}",
    ord(u'\N{LATIN SMALL LETTER N WITH ACUTE}'): "{\\'n}",
    ord(u'\N{LATIN CAPITAL LETTER N WITH CEDILLA}'): '{\\c{N}}',
    ord(u'\N{LATIN SMALL LETTER N WITH CEDILLA}'): '{\\c{n}}',
    ord(u'\N{LATIN CAPITAL LETTER N WITH CARON}'): '{\\v{N}}',
    ord(u'\N{LATIN SMALL LETTER N WITH CARON}'): '{\\v{n}}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH MACRON}'): '{\\=O}',
    ord(u'\N{LATIN SMALL LETTER O WITH MACRON}'): '{\\=o}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH BREVE}'): '{\\u{O}}',
    ord(u'\N{LATIN SMALL LETTER O WITH BREVE}'): '{\\u{o}}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH DOUBLE ACUTE}'): '{\\H{O}}',
    ord(u'\N{LATIN SMALL LETTER O WITH DOUBLE ACUTE}'): '{\\H{o}}',
    ord(u'\N{LATIN CAPITAL LIGATURE OE}'): '{\\OE}',
    ord(u'\N{LATIN SMALL LIGATURE OE}'): '{\\oe}',
    ord(u'\N{LATIN CAPITAL LETTER R WITH ACUTE}'): "{\\'R}",
    ord(u'\N{LATIN SMALL LETTER R WITH ACUTE}'): "{\\'r}",
    ord(u'\N{LATIN CAPITAL LETTER R WITH CEDILLA}'): '{\\c{R}}',
    ord(u'\N{LATIN SMALL LETTER R WITH CEDILLA}'): '{\\c{r}}',
    ord(u'\N{LATIN CAPITAL LETTER R WITH CARON}'): '{\\v{R}}',
    ord(u'\N{LATIN SMALL LETTER R WITH CARON}'): '{\\v{r}}',
    ord(u'\N{LATIN CAPITAL LETTER S WITH ACUTE}'): "{\\'S}",
    ord(u'\N{LATIN SMALL LETTER S WITH ACUTE}'): "{\\'s}",
    ord(u'\N{LATIN CAPITAL LETTER S WITH CIRCUMFLEX}'): '{\\^S}',
    ord(u'\N{LATIN SMALL LETTER S WITH CIRCUMFLEX}'): '{\\^s}',
    ord(u'\N{LATIN CAPITAL LETTER S WITH CEDILLA}'): '{\\c{S}}',
    ord(u'\N{LATIN SMALL LETTER S WITH CEDILLA}'): '{\\c{s}}',
    ord(u'\N{LATIN CAPITAL LETTER S WITH CARON}'): '{\\v{S}}',
    ord(u'\N{LATIN SMALL LETTER S WITH CARON}'): '{\\v{s}}',
    ord(u'\N{LATIN CAPITAL LETTER T WITH CEDILLA}'): '{\\c{T}}',
    ord(u'\N{LATIN SMALL LETTER T WITH CEDILLA}'): '{\\c{t}}',
    ord(u'\N{LATIN CAPITAL LETTER T WITH CARON}'): '{\\v{T}}',
    ord(u'\N{LATIN SMALL LETTER T WITH CARON}'): '{\\v{t}}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH TILDE}'): '{\\~U}',
    ord(u'\N{LATIN SMALL LETTER U WITH TILDE}'): '{\\~u}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH MACRON}'): '{\\=U}',
    ord(u'\N{LATIN SMALL LETTER U WITH MACRON}'): '{\\=u}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH BREVE}'): '{\\u{U}}',
    ord(u'\N{LATIN SMALL LETTER U WITH BREVE}'): '{\\u{u}}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH RING ABOVE}'): '{\\r{U}}',
    ord(u'\N{LATIN SMALL LETTER U WITH RING ABOVE}'): '{\\r{u}}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH DOUBLE ACUTE}'): '{\\H{U}}',
    ord(u'\N{LATIN SMALL LETTER U WITH DOUBLE ACUTE}'): '{\\H{u}}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH OGONEK}'): '{\\c{U}}',
    ord(u'\N{LATIN SMALL LETTER U WITH OGONEK}'): '{\\c{u}}',
    ord(u'\N{LATIN CAPITAL LETTER W WITH CIRCUMFLEX}'): '{\\^W}',
    ord(u'\N{LATIN SMALL LETTER W WITH CIRCUMFLEX}'): '{\\^w}',
    ord(u'\N{LATIN CAPITAL LETTER Y WITH CIRCUMFLEX}'): '{\\^Y}',
    ord(u'\N{LATIN SMALL LETTER Y WITH CIRCUMFLEX}'): '{\\^y}',
    ord(u'\N{LATIN CAPITAL LETTER Y WITH DIAERESIS}'): '{\\"Y}',
    ord(u'\N{LATIN CAPITAL LETTER Z WITH ACUTE}'): "{\\'Z}",
    ord(u'\N{LATIN SMALL LETTER Z WITH ACUTE}'): "{\\'Z}",
    ord(u'\N{LATIN CAPITAL LETTER Z WITH DOT ABOVE}'): '{\\.Z}',
    ord(u'\N{LATIN SMALL LETTER Z WITH DOT ABOVE}'): '{\\.Z}',
    ord(u'\N{LATIN CAPITAL LETTER Z WITH CARON}'): '{\\v{Z}}',
    ord(u'\N{LATIN SMALL LETTER Z WITH CARON}'): '{\\v{z}}',
    ord(u'\N{LATIN CAPITAL LETTER DZ WITH CARON}'): '{D\\v{Z}}',
    ord(u'\N{LATIN CAPITAL LETTER D WITH SMALL LETTER Z WITH CARON}'): '{D\\v{z}}',
    ord(u'\N{LATIN SMALL LETTER DZ WITH CARON}'): '{d\\v{z}}',
    ord(u'\N{LATIN CAPITAL LETTER LJ}'): '{LJ}',
    ord(u'\N{LATIN CAPITAL LETTER L WITH SMALL LETTER J}'): '{Lj}',
    ord(u'\N{LATIN SMALL LETTER LJ}'): '{lj}',
    ord(u'\N{LATIN CAPITAL LETTER NJ}'): '{NJ}',
    ord(u'\N{LATIN CAPITAL LETTER N WITH SMALL LETTER J}'): '{Nj}',
    ord(u'\N{LATIN SMALL LETTER NJ}'): '{nj}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH CARON}'): '{\\v{A}}',
    ord(u'\N{LATIN SMALL LETTER A WITH CARON}'): '{\\v{a}}',
    ord(u'\N{LATIN CAPITAL LETTER I WITH CARON}'): '{\\v{I}}',
    ord(u'\N{LATIN SMALL LETTER I WITH CARON}'): '{\\v\\i}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH CARON}'): '{\\v{O}}',
    ord(u'\N{LATIN SMALL LETTER O WITH CARON}'): '{\\v{o}}',
    ord(u'\N{LATIN CAPITAL LETTER U WITH CARON}'): '{\\v{U}}',
    ord(u'\N{LATIN SMALL LETTER U WITH CARON}'): '{\\v{u}}',
    ord(u'\N{LATIN CAPITAL LETTER G WITH CARON}'): '{\\v{G}}',
    ord(u'\N{LATIN SMALL LETTER G WITH CARON}'): '{\\v{g}}',
    ord(u'\N{LATIN CAPITAL LETTER K WITH CARON}'): '{\\v{K}}',
    ord(u'\N{LATIN SMALL LETTER K WITH CARON}'): '{\\v{k}}',
    ord(u'\N{LATIN CAPITAL LETTER O WITH OGONEK}'): '{\\c{O}}',
    ord(u'\N{LATIN SMALL LETTER O WITH OGONEK}'): '{\\c{o}}',
    ord(u'\N{LATIN SMALL LETTER J WITH CARON}'): '{\\v\\j}',
    ord(u'\N{LATIN CAPITAL LETTER DZ}'): '{DZ}',
    ord(u'\N{LATIN CAPITAL LETTER D WITH SMALL LETTER Z}'): '{Dz}',
    ord(u'\N{LATIN SMALL LETTER DZ}'): '{dz}',
    ord(u'\N{LATIN CAPITAL LETTER G WITH ACUTE}'): "{\\'G}",
    ord(u'\N{LATIN SMALL LETTER G WITH ACUTE}'): "{\\'g}",
    ord(u'\N{LATIN CAPITAL LETTER AE WITH ACUTE}'): "{\\'\\AE}",
    ord(u'\N{LATIN SMALL LETTER AE WITH ACUTE}'): "{\\'\\ae}",
    ord(u'\N{LATIN CAPITAL LETTER O WITH STROKE AND ACUTE}'): "{\\'\\O}",
    ord(u'\N{LATIN SMALL LETTER O WITH STROKE AND ACUTE}'): "{\\'\\o}",
    ord(u'\N{PARTIAL DIFFERENTIAL}'): '{\\mbox{$\\partial$}}',
    ord(u'\N{N-ARY PRODUCT}'): '{\\mbox{$\\prod$}}',
    ord(u'\N{N-ARY SUMMATION}'): '{\\mbox{$\\sum$}}',
    ord(u'\N{SQUARE ROOT}'): '{\\mbox{$\\surd$}}',
    ord(u'\N{INFINITY}'): '{\\mbox{$\\infty$}}',
    ord(u'\N{INTEGRAL}'): '{\\mbox{$\\int$}}',
    ord(u'\N{INTERSECTION}'): '{\\mbox{$\\cap$}}',
    ord(u'\N{UNION}'): '{\\mbox{$\\cup$}}',
    ord(u'\N{RIGHTWARDS ARROW}'): '{\\mbox{$\\rightarrow$}}',
    ord(u'\N{RIGHTWARDS DOUBLE ARROW}'): '{\\mbox{$\\Rightarrow$}}',
    ord(u'\N{LEFTWARDS ARROW}'): '{\\mbox{$\\leftarrow$}}',
    ord(u'\N{LEFTWARDS DOUBLE ARROW}'): '{\\mbox{$\\Leftarrow$}}',
    ord(u'\N{LOGICAL OR}'): '{\\mbox{$\\vee$}}',
    ord(u'\N{LOGICAL AND}'): '{\\mbox{$\\wedge$}}',
    ord(u'\N{ALMOST EQUAL TO}'): '{\\mbox{$\\approx$}}',
    ord(u'\N{NOT EQUAL TO}'): '{\\mbox{$\\neq$}}',
    ord(u'\N{LESS-THAN OR EQUAL TO}'): '{\\mbox{$\\leq$}}',
    ord(u'\N{GREATER-THAN OR EQUAL TO}'): '{\\mbox{$\\geq$}}',
    ord(u'\N{MODIFIER LETTER CIRCUMFLEX ACCENT}'): '{\\^{}}',
    ord(u'\N{CARON}'): '{\\v{}}',
    ord(u'\N{BREVE}'): '{\\u{}}',
    ord(u'\N{DOT ABOVE}'): '{\\.{}}',
    ord(u'\N{RING ABOVE}'): '{\\r{}}',
    ord(u'\N{OGONEK}'): '{\\c{}}',
    ord(u'\N{SMALL TILDE}'): '{\\~{}}',
    ord(u'\N{DOUBLE ACUTE ACCENT}'): '{\\H{}}',
    ord(u'\N{LATIN SMALL LIGATURE FI}'): '{fi}',
    ord(u'\N{LATIN SMALL LIGATURE FL}'): '{fl}',
    ord(u'\N{GREEK SMALL LETTER ALPHA}'): '{\\mbox{$\\alpha$}}',
    ord(u'\N{GREEK SMALL LETTER BETA}'): '{\\mbox{$\\beta$}}',
    ord(u'\N{GREEK SMALL LETTER GAMMA}'): '{\\mbox{$\\gamma$}}',
    ord(u'\N{GREEK SMALL LETTER DELTA}'): '{\\mbox{$\\delta$}}',
    ord(u'\N{GREEK SMALL LETTER EPSILON}'): '{\\mbox{$\\epsilon$}}',
    ord(u'\N{GREEK SMALL LETTER ZETA}'): '{\\mbox{$\\zeta$}}',
    ord(u'\N{GREEK SMALL LETTER ETA}'): '{\\mbox{$\\eta$}}',
    ord(u'\N{GREEK SMALL LETTER THETA}'): '{\\mbox{$\\theta$}}',
    ord(u'\N{GREEK SMALL LETTER IOTA}'): '{\\mbox{$\\iota$}}',
    ord(u'\N{GREEK SMALL LETTER KAPPA}'): '{\\mbox{$\\kappa$}}',
    ord(u'\N{GREEK SMALL LETTER LAMDA}'): '{\\mbox{$\\lambda$}}', # NO B??!?
    ord(u'\N{GREEK SMALL LETTER MU}'): '{\\mbox{$\\mu$}}',
    ord(u'\N{GREEK SMALL LETTER NU}'): '{\\mbox{$\\nu$}}',
    ord(u'\N{GREEK SMALL LETTER XI}'): '{\\mbox{$\\xi$}}',
    ord(u'\N{GREEK SMALL LETTER OMICRON}'): '{\\mbox{$\\omicron$}}',
    ord(u'\N{GREEK SMALL LETTER PI}'): '{\\mbox{$\\pi$}}',
    ord(u'\N{GREEK SMALL LETTER RHO}'): '{\\mbox{$\\rho$}}',
    ord(u'\N{GREEK SMALL LETTER SIGMA}'): '{\\mbox{$\\sigma$}}',
    ord(u'\N{GREEK SMALL LETTER TAU}'): '{\\mbox{$\\tau$}}',
    ord(u'\N{GREEK SMALL LETTER UPSILON}'): '{\\mbox{$\\upsilon$}}',
    ord(u'\N{GREEK SMALL LETTER PHI}'): '{\\mbox{$\\phi$}}',
    ord(u'\N{GREEK SMALL LETTER CHI}'): '{\\mbox{$\\chi$}}',
    ord(u'\N{GREEK SMALL LETTER PSI}'): '{\\mbox{$\\psi$}}',
    ord(u'\N{GREEK SMALL LETTER OMEGA}'): '{\\mbox{$\\omega$}}',
    ord(u'\N{GREEK CAPITAL LETTER ALPHA}'): '{\\mbox{$\\Alpha$}}',
    ord(u'\N{GREEK CAPITAL LETTER BETA}'): '{\\mbox{$\\Beta$}}',
    ord(u'\N{GREEK CAPITAL LETTER GAMMA}'): '{\\mbox{$\\Gamma$}}',
    ord(u'\N{GREEK CAPITAL LETTER DELTA}'): '{\\mbox{$\\Delta$}}',
    ord(u'\N{GREEK CAPITAL LETTER EPSILON}'): '{\\mbox{$\\Epsilon$}}',
    ord(u'\N{GREEK CAPITAL LETTER ZETA}'): '{\\mbox{$\\Zeta$}}',
    ord(u'\N{GREEK CAPITAL LETTER ETA}'): '{\\mbox{$\\Eta$}}',
    ord(u'\N{GREEK CAPITAL LETTER THETA}'): '{\\mbox{$\\Theta$}}',
    ord(u'\N{GREEK CAPITAL LETTER IOTA}'): '{\\mbox{$\\Iota$}}',
    ord(u'\N{GREEK CAPITAL LETTER KAPPA}'): '{\\mbox{$\\Kappa$}}',
    ord(u'\N{GREEK CAPITAL LETTER LAMDA}'): '{\\mbox{$\\Lambda$}}', # DITTO
    ord(u'\N{GREEK CAPITAL LETTER MU}'): '{\\mbox{$\\Mu$}}',
    ord(u'\N{GREEK CAPITAL LETTER NU}'): '{\\mbox{$\\Nu$}}',
    ord(u'\N{GREEK CAPITAL LETTER XI}'): '{\\mbox{$\\Xi$}}',
    ord(u'\N{GREEK CAPITAL LETTER OMICRON}'): '{\\mbox{$\\Omicron$}}',
    ord(u'\N{GREEK CAPITAL LETTER PI}'): '{\\mbox{$\\Pi$}}',
    ord(u'\N{GREEK CAPITAL LETTER RHO}'): '{\\mbox{$\\Rho$}}',
    ord(u'\N{GREEK CAPITAL LETTER SIGMA}'): '{\\mbox{$\\Sigma$}}',
    ord(u'\N{GREEK CAPITAL LETTER TAU}'): '{\\mbox{$\\Tau$}}',
    ord(u'\N{GREEK CAPITAL LETTER UPSILON}'): '{\\mbox{$\\Upsilon$}}',
    ord(u'\N{GREEK CAPITAL LETTER PHI}'): '{\\mbox{$\\Phi$}}',
    ord(u'\N{GREEK CAPITAL LETTER CHI}'): '{\\mbox{$\\Chi$}}',
    ord(u'\N{GREEK CAPITAL LETTER PSI}'): '{\\mbox{$\\Psi$}}',
    ord(u'\N{GREEK CAPITAL LETTER OMEGA}'): '{\\mbox{$\\Omega$}}',
    ord(u'\N{COPYRIGHT SIGN}'): '{\\copyright}',
    ord(u'\N{LATIN CAPITAL LETTER A WITH ACUTE}'): "{\\'A}",
    ord(u'\N{LATIN CAPITAL LETTER I WITH ACUTE}'): "{\\'I}",
    ord(u'\N{HORIZONTAL ELLIPSIS}'): '{\\ldots}',
    ord(u'\N{TRADE MARK SIGN}'): '{\\mbox{$^\\mbox{TM}$}}',
}
for _i in range(0x0020):
    if _i not in latex_equivalents:
        latex_equivalents[_i] = ''
for _i in range(0x0020,0x007f):
    if _i not in latex_equivalents:
        latex_equivalents[_i] = chr(_i)

# Characters that should be ignored and not output in tokenization
_ignore = set([chr(i) for i in range(32)+[127]]) - set('\t\n\r')

# Regexp of chars not in blacklist, for quick start of tokenize
_stoppers = re.compile('[\x00-\x1f!$\\-?\\{~\\\\`\']')

_blacklist = set(' \n\r')
_blacklist.add(None)    # shortcut candidate generation at end of data

# Construction of inverse translation table
_l2u = {
    '\ ':ord(' ')   # unexpanding space makes no sense in non-TeX contexts
}

for _tex in latex_equivalents:
    if _tex <= 0x0020 or (_tex <= 0x007f and len(latex_equivalents[_tex]) <= 1):
        continue    # boring entry
    _toks = tuple(_tokenize(latex_equivalents[_tex]))
    if _toks[0] == '{' and _toks[-1] == '}':
        _toks = _toks[1:-1]
    if _toks[0].isalpha():
        continue    # don't turn ligatures into single chars
    if len(_toks) == 1 and (_toks[0] == "'" or _toks[0] == "`"):
        continue    # don't turn ascii quotes into curly quotes
    if _toks[0] == '\\mbox' and _toks[1] == '{' and _toks[-1] == '}':
        _toks = _toks[2:-1]
    if len(_toks) == 4 and _toks[1] == '{' and _toks[3] == '}':
        _toks = (_toks[0],_toks[2])
    if len(_toks) == 1:
        _toks = _toks[0]
    _l2u[_toks] = _tex

# Shortcut candidate generation for certain useless candidates:
# a character is in _blacklist if it can not be at the start
# of any translation in _l2u.  We use this to quickly skip through
# such characters before getting to more difficult-translate parts.
# _blacklist is defined several lines up from here because it must
# be defined in order to call _tokenize, however it is safe to
# delay filling it out until now.

for i in range(0x0020,0x007f):
    _blacklist.add(chr(i))
_blacklist.remove('{')
_blacklist.remove('$')
for candidate in _l2u:
    if isinstance(candidate,tuple):
        if not candidate or not candidate[0]:
            continue
        firstchar = candidate[0][0]
    else:
        firstchar = candidate[0]
    _blacklist.discard(firstchar)
