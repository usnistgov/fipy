#
# epydoc package file
#
# A python documentation Module
# Edward Loper
#
# $Id$
#

"""
Markup language support for docstrings.  Each submodule defines a
parser for a single markup language.  These parsers convert an
object's docstring to a L{ParsedDocstring}, a standard intermediate
representation that can be used to generate output.
C{ParsedDocstring}s support the following operations:
  - output generation (L{to_plaintext()<ParsedDocstring.to_plaintext>},
    L{to_html()<ParsedDocstring.to_html>}, and
    L{to_latex()<ParsedDocstring.to_latex>}).
  - Summarization (L{summary()<ParsedDocstring.summary>}).
  - Field extraction (L{split_fields()<ParsedDocstring.split_fields>}).
  - Index term extraction (L{index_terms()<ParsedDocstring.index_terms>}.

The L{parse()} function provides a single interface to the
C{epydoc.markup} package: it takes a docstring and the name of a
markup language; delegates to the appropriate parser; and returns the
parsed docstring (along with any errors or warnings that were
generated).

The C{ParsedDocstring} output generation methods (C{to_M{format}()})
use a L{DocstringLinker} to link the docstring output with the rest of
the documentation that epydoc generates.  C{DocstringLinker}s are
currently responsible for translating two kinds of crossreference:
  - index terms (L{translate_indexterm()
    <DocstringLinker.translate_indexterm>}).
  - identifier crossreferences (L{translate_identifier_xref()
    <DocstringLinker.translate_identifier_xref>}).

A parsed docstring's fields can be extracted using the
L{ParsedDocstring.split_fields()} method.  This method divides a
docstring into its main body and a list of L{Field}s, each of which
encodes a single field.  The field's bodies are encoded as
C{ParsedDocstring}s.

Markup errors are represented using L{ParseError}s.  These exception
classes record information about the cause, location, and severity of
each error.

The C{epydoc.markup} module also defines several utility functions,
such as L{wordwrap}, L{plaintext_to_latex}, and L{plaintext_to_html},
which are used by several different markup language parsers.

@sort: parse, ParsedDocstring, Field, DocstringLinker
@group Errors and Warnings: ParseError
@group Utility Functions: wordwrap, plaintet_to_html, plaintext_to_latex,
    parse_type_of
@var SCRWIDTH: The default width with which text will be wrapped
      when formatting the output of the parser.
@type SCRWIDTH: C{int}
@var _parse_warnings: Used by L{_parse_warn}.
"""
__docformat__ = 'epytext en'

import re, types, sys
from epydoc.imports import import_module

##################################################
## Contents
##################################################
#
# 1. parse() dispatcher
# 2. ParsedDocstring abstract base class
# 3. Field class
# 4. Docstring Linker
# 5. ParseError exceptions
# 6. Misc helpers
#

##################################################
## Dispatcher
##################################################

def parse(docstring, markup='plaintext', errors=None, **options):
    """
    Parse the given docstring, and use it to construct a
    C{ParsedDocstring}.  If any fatal C{ParseError}s are encountered
    while parsing the docstring, then the docstring will be rendered
    as plaintext, instead.

    @type docstring: C{string}
    @param docstring: The docstring to encode.
    @type markup: C{string}
    @param markup: The name of the markup language that is used by
        the docstring.  If the markup language is not supported, then
        the docstring will be treated as plaintext.  The markup name
        is case-insensitive.
    @param errors: A list where any errors generated during parsing
        will be stored.  If no list is specified, then fatal errors
        will generate exceptions, and non-fatal errors will be
        ignored.
    @type errors: C{list} of L{ParseError}
    @rtype: L{ParsedDocstring}
    @return: A L{ParsedDocstring} that encodes the contents of
        C{docstring}.
    @raise ParseError: If C{errors} is C{None} and an error is
        encountered while parsing.
    """
    # Initialize errors list.
    raise_on_error = (errors is None)
    if errors == None: errors = []

    # Normalize the markup language name.
    markup = markup.lower()

    # Is the markup language valid?
    if not re.match(r'\w+', markup):
        _parse_warn('Warning: Bad markup language name: %s' % markup)
        return plaintext.parse_docstring(docstring, errors, **options)

    # Is the markup language supported?
    try: exec('from epydoc.markup.%s import parse_docstring' % markup)
    except:
        _parse_warn('Warning: Unsupported markup language: %s' % markup)
        return plaintext.parse_docstring(docstring, errors, **options)

    # Parse the docstring.
    try: parsed_docstring = parse_docstring(docstring, errors, **options)
    except KeyboardInterrupt: raise
    except Exception, e:
        errors.append(ParseError('Internal error: %s' % e))
        return plaintext.parse_docstring(docstring, errors, **options)

    # Check for fatal errors.
    fatal_errors = [e for e in errors if e.is_fatal()]
    if fatal_errors and raise_on_error: raise fatal_errors[0]
    if fatal_errors:
        return plaintext.parse_docstring(docstring, errors, **options)

    return parsed_docstring

# only issue each warning once:
_parse_warnings = {}
def _parse_warn(estr):
    """
    Print a warning message.  If the given error has already been
    printed, then do nothing.
    """
    global _parse_warnings
    if _parse_warnings.has_key(estr): return
    _parse_warnings[estr] = 1
    if sys.stderr.softspace: print >>sys.stderr
    print >>sys.stderr, estr

##################################################
## ParsedDocstring
##################################################
class ParsedDocstring:
    """
    A standard intermediate representation for parsed docstrings that
    can be used to generate output.  Parsed docstrings are produced by
    markup parsers (such as L{epytext.parse} or L{javadoc.parse}).
    C{ParsedDocstring}s support several kinds of operation:    
      - output generation (L{to_plaintext()}, L{to_html()}, and
        L{to_latex()}).
      - Summarization (L{summary()}).
      - Field extraction (L{split_fields()}).
      - Index term extraction (L{index_terms()}.

    The output generation methods (C{to_M{format}()}) use a
    L{DocstringLinker} to link the docstring output with the rest
    of the documentation that epydoc generates.

    Subclassing
    ===========
    The only method that a subclass is I{required} to implement is
    L{to_plaintext()}; but it is often useful to override the other
    methods.  The default behavior of each method is described below:
      - C{to_I{format}}: Calls C{to_plaintext}, and uses the string it
        returns to generate verbatim output.
      - C{summary}: Returns C{self} (i.e., the entire docstring).
      - C{split_fields}: Returns C{(self, [])} (i.e., extracts no
        fields).
      - C{index_terms}: Returns C{[]} (i.e., extracts no index terms).

    If and when epydoc adds more output formats, new C{to_I{format}}
    methods will be added to this base class; but they will always
    be given a default implementation.
    """
    def split_fields(self, errors=None):
        """
        Split this docstring into its body and its fields.
        
        @return: A tuple C{(M{body}, M{fields})}, where C{M{body}} is
            the main body of this docstring, and C{M{fields}} is a list
            of its fields.
        @rtype: C{(L{ParsedDocstring}, list of L{Field})}
        @param errors: A list where any errors generated during
            splitting will be stored.  If no list is specified, then
            errors will be ignored.
        @type errors: C{list} of L{ParseError}
        """
        # Default behavior:
        return self, []

    def summary(self):
        """
        @return: A short summary of this docstring.  Typically, the
            summary consists of the first sentence of the docstring.
        @rtype: L{ParsedDocstring}
        """
        # Default behavior:
        return self

    def concatenate(self, other):
        """
        @return: A new parsed docstring containing the concatination
            of this docstring and C{other}.
        @raise ValueError: If the two parsed docstrings are
            incompatible.
        """
        # Default behavior:
        raise ValueError, 'Could not concatenate docstrings'

    def __add__(self, other): return self.concatenate(other)

    def to_html(self, docstring_linker, **options):
        """
        Translate this docstring to HTML.

        @param docstring_linker: An HTML translator for crossreference
            links into and out of the docstring.
        @type docstring_linker: L{DocstringLinker}
        @param options: Any extra options for the output.  Unknown
            options are ignored.
        @return: An HTML fragment that encodes this docstring.
        @rtype: C{string}
        """
        # Default behavior:
        plaintext = plaintext_to_html(self.to_plaintext(docstring_linker))
        return '<pre class="literalblock">\n%s\n</pre>\n' % plaintext

    def to_latex(self, docstring_linker, **options):
        """
        Translate this docstring to LaTeX.
        
        @param docstring_linker: A LaTeX translator for crossreference
            links into and out of the docstring.
        @type docstring_linker: L{DocstringLinker}
        @param options: Any extra options for the output.  Unknown
            options are ignored.
        @return: A LaTeX fragment that encodes this docstring.
        @rtype: C{string}
        """
        # Default behavior:
        plaintext = plaintext_to_latex(self.to_plaintext(docstring_linker))
        return '\\begin{alltt}\n%s\\end{alltt}\n\n' % plaintext

    def to_plaintext(self, docstring_linker, **options):
        """
        Translate this docstring to plaintext.
        
        @param docstring_linker: A plaintext translator for
            crossreference links into and out of the docstring.
        @type docstring_linker: L{DocstringLinker}
        @param options: Any extra options for the output.  Unknown
            options are ignored.
        @return: A plaintext fragment that encodes this docstring.
        @rtype: C{string}
        """
        raise NotImplementedError, 'ParsedDocstring.to_plaintext()'

    def index_terms(self):
        """
        @return: The list of index terms that are defined in this
            docstring.  Each of these items will be added to the index
            page of the documentation.
        @rtype: C{list} of C{ParsedDocstring}
        """
        # Default behavior:
        return []

##################################################
## Fields
##################################################
class Field:
    """
    The contents of a docstring's field.  Docstring fields are used
    to describe specific aspects of an object, such as a parameter of
    a function or the author of a module.  Each field consists of a
    tag, an optional argument, and a body:
      - The tag specifies the type of information that the field
        encodes.
      - The argument specifies the object that the field describes.
        The argument may be C{None} or a C{string}.
      - The body contains the field's information.

    Tags are automatically downcased and stripped; and arguments are
    automatically stripped.
    """
    def __init__(self, tag, arg, body):
        self._tag = tag.lower().strip()
        if arg is None: self._arg = None
        else: self._arg = arg.strip()
        self._body = body

    def tag(self):
        """
        @return: This field's tag.
        @rtype: C{string}
        """
        return self._tag

    def arg(self):
        """
        @return: This field's argument, or C{None} if this field has
            no argument.
        @rtype: C{string} or C{None}
        """
        return self._arg

    def body(self):
        """
        @return: This field's body.
        @rtype: L{ParsedDocstring}
        """
        return self._body

    def __repr__(self):
        if self._arg is None:
            return '<Field @%s: ...>' % self._tag
        else:
            return '<Field @%s %s: ...>' % (self._tag, self._arg)

##################################################
## Docstring Linker (resolves crossreferences)
##################################################
class DocstringLinker: 
    """
    A translator for crossreference links into and out of a
    C{ParsedDocstring}.  C{DocstringLinker} is used by
    C{ParsedDocstring} to convert these crossreference links into
    appropriate output formats.  For example,
    C{DocstringLinker.to_html} expects a C{DocstringLinker} that
    converts crossreference links to HTML.
    """
    def translate_indexterm(self, indexterm):
        """
        Translate an index term to the appropriate output format.  The
        output will typically include a crossreference anchor.

        @type indexterm: L{ParsedDocstring}
        @param indexterm: The index term to translate.
        @rtype: C{string}
        @return: The translated index term.
        """
        raise NotImplementedError, 'DocstringLinker.translate_indexterm()'

    def translate_identifier_xref(self, identifier, label=None):
        """
        Translate a crossreference link to a Python identifier to the
        appropriate output format.  The output will typically include
        a reference or pointer to the crossreference target.

        @type identifier: C{string}
        @param identifier: The name of the Python identifier that
            should be linked to.
        @type label: C{string} or C{None}
        @param label: The label that should be used for the identifier,
            if it's different from the name of the identifier.
        @rtype: C{string}
        @return: The translated crossreference link.
        """
        raise NotImplementedError, 'DocstringLinker.translate_xref()'

##################################################
## ParseError exceptions
##################################################

class ParseError(Exception):
    """
    The base class for errors generated while parsing docstrings.

    @ivar _linenum: The line on which the error occured within the
        docstring.  The linenum of the first line is 0.
    @type _linenum: C{int}
    @ivar _offset: The line number where the docstring begins.  This
        offset is added to C{_linenum} when displaying the line number
        of the error.  Default value: 1.
    @type _offset: C{int}
    @ivar _descr: A description of the error.
    @type _descr: C{string}
    @ivar _fatal: True if this is a fatal error.
    @type _fatal: C{boolean}
    """
    def __init__(self, descr, linenum=None, is_fatal=1):
        """
        @type descr: C{string}
        @param descr: A description of the error.
        @type linenum: C{int}
        @param linenum: The line on which the error occured within
            the docstring.  The linenum of the first line is 0.
        @type is_fatal: C{boolean}
        @param is_fatal: True if this is a fatal error.
        """
        self._descr = descr
        self._linenum = linenum
        self._fatal = is_fatal
        self._offset = 1
                 
    def is_fatal(self):
        """
        @return: true if this is a fatal error.  If an error is fatal,
            then epydoc should ignore the output of the parser, and
            parse the docstring as plaintext.
        @rtype: C{boolean}
        """
        return self._fatal

    def linenum(self):
        """
        @return: The line number on which the error occured (including
        any offset).  If the line number is unknown, then return
        C{None}.
        @rtype: C{int} or C{None}
        """
        if self._linenum is None: return None
        else: return self._offset + self._linenum

    def set_linenum_offset(self, offset):
        """
        Set the line number offset for this error.  This offset is the
        line number where the docstring begins.  This offset is added
        to C{_linenum} when displaying the line number of the error.

        @param offset: The new line number offset.
        @type offset: C{int}
        @rtype: C{None}
        """
        self._offset = offset
    
    def __str__(self):
        """
        Return a string representation of this C{ParseError}.  This
        multi-line string contains a description of the error, and
        specifies where it occured.
        
        @return: the informal representation of this C{ParseError}.
        @rtype: C{string}
        """
        if self._linenum is not None:
            str = '%5s: ' % ('L'+`self._linenum+self._offset`)
        else:
            str = '     - '
        if self._fatal:
            str += 'Error: '
        else:
            str += 'Warning: '
        return str + wordwrap(self._descr, 7, startindex=len(str))[:-1]
    
    def __repr__(self):
        """
        Return the formal representation of this C{ParseError}.
        C{ParseError}s have formal representations of the form::
           <ParseError on line 12>

        @return: the formal representation of this C{ParseError}.
        @rtype: C{string}
        """
        if self._linenum is None:
            return '<ParseError on line %d' % self._offset
        else:
            return '<ParseError on line %d>' % (self._linenum+self._offset)

    def __cmp__(self, other):
        """
        Compare two C{ParseError}s, based on their line number.
          - Return -1 if C{self.linenum<other.linenum}
          - Return +1 if C{self.linenum>other.linenum}
          - Return 0 if C{self.linenum==other.linenum}.
        The return value is undefined if C{other} is not a
        ParseError.

        @rtype: C{int}
        """
        if not isinstance(other, ParseError): return -1000
        return cmp(self._linenum+self._offset,
                   other._linenum+other._offset)

##################################################
## Misc helpers
##################################################
# These are used by multiple markup parsers

# Default screen width, for word-wrapping
SCRWIDTH = 73

def wordwrap(str, indent=0, right=SCRWIDTH, startindex=0):
    """
    Word-wrap the given string.  All sequences of whitespace are
    converted into spaces, and the string is broken up into lines,
    where each line begins with C{indent} spaces, followed by one or
    more (space-deliniated) words whose length is less than
    C{right-indent}.  If a word is longer than C{right-indent}
    characters, then it is put on its own line.

    @param str: The string that should be word-wrapped.
    @type str: C{int}
    @param indent: The left margin of the string.  C{indent} spaces
        will be inserted at the beginning of every line.
    @type indent: C{int}
    @param right: The right margin of the string.
    @type right: C{int}
    @type startindex: C{int}
    @param startindex: The index at which the first line starts.  This
        is useful if you want to include other contents on the first
        line. 
    @return: A word-wrapped version of C{str}.
    @rtype: C{string}
    """
    words = str.split()
    out_str = ' '*(indent-startindex)
    charindex = max(indent, startindex)
    for word in words:
        if charindex+len(word) > right and charindex > 0:
            out_str += '\n' + ' '*indent
            charindex = indent
        out_str += word+' '
        charindex += len(word)+1
    return out_str.rstrip()+'\n'

def plaintext_to_html(str):
    """
    @return: An HTML string that encodes the given plaintext string.
    In particular, special characters (such as C{'<'} and C{'&'})
    are escaped.
    @rtype: C{string}
    """
    str = str.replace('&', '&amp;').replace('"', '&quot;')
    str = str.replace('<', '&lt;').replace('>', '&gt;')
    return str.replace('@', '&#64;')

def plaintext_to_latex(str, nbsp=0, breakany=0):
    """
    @return: A LaTeX string that encodes the given plaintext string.
    In particular, special characters (such as C{'$'} and C{'_'})
    are escaped, and tabs are expanded.
    @rtype: C{string}
    @param breakany: Insert hyphenation marks, so that LaTeX can
    break the resulting string at any point.  This is useful for
    small boxes (e.g., the type box in the variable list table).
    @param nbsp: Replace every space with a non-breaking space
    (C{'~'}).
    """
    # These get converted to hyphenation points later
    if breakany: str = re.sub('(.)', '\\1\1', str)

    # These get converted to \textbackslash later.
    str = str.replace('\\', '\0')

    # Expand tabs
    str = str.expandtabs()

    # These elements need to be backslashed.
    str = re.sub(r'([#$&%_\${}])', r'\\\1', str)

    # These elements have special names.
    str = str.replace('|', '{\\textbar}')
    str = str.replace('<', '{\\textless}')
    str = str.replace('>', '{\\textgreater}')
    str = str.replace('^', '{\\textasciicircum}')
    str = str.replace('~', '{\\textasciitilde}')
    str = str.replace('\0', r'{\textbackslash}')

    # replace spaces with non-breaking spaces
    if nbsp: str = str.replace(' ', '~')

    # Convert \1's to hyphenation points.
    if breakany: str = str.replace('\1', r'\-')
    
    return str

def parse_type_of(obj):
    """
    @return: A C{ParsedDocstring} that encodes the type of the given
    object.
    @rtype: L{ParsedDocstring}
    @param obj: The object whose type should be returned as DOM document.
    @type obj: any
    """
    # This is a bit hackish; oh well. :)
    from epydoc.markup.epytext import ParsedEpytextDocstring
    from xml.dom.minidom import Document
    doc = Document()
    epytext = doc.createElement('epytext')
    para = doc.createElement('para')
    doc.appendChild(epytext)
    epytext.appendChild(para)
    
    if type(obj) is types.InstanceType:
        link = doc.createElement('link')
        name = doc.createElement('name')
        target = doc.createElement('target')
        para.appendChild(link)
        link.appendChild(name)
        link.appendChild(target)
        name.appendChild(doc.createTextNode(str(obj.__class__.__name__)))
        target.appendChild(doc.createTextNode(str(obj.__class__)))        
    else:
        code = doc.createElement('code')
        para.appendChild(code)
        code.appendChild(doc.createTextNode(type(obj).__name__))
    return ParsedEpytextDocstring(doc)

##################################################
## Sub-module Imports
##################################################
# By default, just import plaintext.  That way we don't have to wait
# for other modules (esp restructuredtext) to load if we're not going
# to use them.

import plaintext
