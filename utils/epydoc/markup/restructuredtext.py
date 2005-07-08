#
# rst.py: ReStructuredText docstring parsing
# Edward Loper
#
# Created [06/28/03 02:52 AM]
# $Id$
#

"""
Epydoc parser for ReStructuredText strings.  ReStructuredText is the
standard markup language used by the Docutils project.
L{parse_docstring()} provides the primary interface to this module; it
returns a L{ParsedRstDocstring}, which supports all of the methods
defined by L{ParsedDocstring}.

L{ParsedRstDocstring} is basically just a L{ParsedDocstring} wrapper
for the L{docutils.nodes.document} class.

Creating C{ParsedRstDocstring}s
===============================
C{ParsedRstDocstring}s are created by the C{parse_document} function,
using the L{docutils.core.publish_string()} method, with the following
helpers:

  - An L{_EpydocReader} is used to capture all error messages as it
    parses the docstring.
  - A L{_DocumentPseudoWriter} is used to extract the document itself,
    without actually writing any output.  The document is saved for
    further processing.  The settings for the writer are copied from
    L{docutils.writers.html4css1.Writer}, since those settings will
    be used when we actually write the docstring to html.

Using C{ParsedRstDocstring}s
============================

C{ParsedRstDocstring}s support all of the methods defined by
C{ParsedDocstring}; but only the following four methods have
non-default behavior:

  - L{to_html()<ParsedRstDocstring.to_html>} uses an
    L{_EpydocHTMLTranslator} to translate the C{ParsedRstDocstring}'s
    document into an HTML segment.
  - L{split_fields()<ParsedRstDocstring.split_fields>} uses a
    L{_SplitFieldsTranslator} to divide the C{ParsedRstDocstring}'s
    document into its main body and its fields.  Special handling
    is done to account for consolidated fields.
  - L{summary()<ParsedRstDocstring.summary>} uses a
    L{_SummaryExtractor} to extract the first sentence from
    the C{ParsedRstDocstring}'s document.
  - L{to_plaintext()<ParsedRstDocstring.to_plaintext>} uses
    C{document.astext()} to convert the C{ParsedRstDocstring}'s
    document to plaintext.

@todo: Add ParsedRstDocstring.to_latex()
@var CONSOLIDATED_FIELDS: A dictionary encoding the set of
'consolidated fields' that can be used.  Each consolidated field is
marked by a single tag, and contains a single bulleted list, where
each list item starts with an identifier, marked as interpreted text
(C{`...`}).  This module automatically splits these consolidated
fields into individual fields.  The keys of C{CONSOLIDATED_FIELDS} are
the names of possible consolidated fields; and the values are the
names of the field tags that should be used for individual entries in
the list.
"""
__docformat__ = 'epytext en'

# Imports
import re
from xml.dom.minidom import *

from docutils.core import publish_string
from docutils.writers import Writer
from docutils.writers.html4css1 import HTMLTranslator, Writer as HTMLWriter
from docutils.writers.latex2e import LaTeXTranslator, Writer as LaTeXWriter
from docutils.readers.standalone import Reader as StandaloneReader
from docutils.utils import new_document
from docutils.nodes import NodeVisitor, Text, SkipChildren
from docutils.nodes import SkipNode, TreeCopyVisitor
from docutils.frontend import OptionParser
import docutils.nodes

from epydoc.markup import *

CONSOLIDATED_FIELDS = {
    'parameters': 'param',
    'arguments': 'arg',
    'exceptions': 'except',
    'variables': 'var',
    'ivariables': 'ivar',
    'cvariables': 'cvar',
    'groups': 'group',
    'types': 'type',
    'keywords': 'keyword',
    }

def parse_docstring(docstring, errors, **options):
    """
    Parse the given docstring, which is formatted using
    ReStructuredText; and return a L{ParsedDocstring} representation
    of its contents.
    @param docstring: The docstring to parse
    @type docstring: C{string}
    @param errors: A list where any errors generated during parsing
        will be stored.
    @type errors: C{list} of L{ParseError}
    @param options: Extra options.  Unknown options are ignored.
        Currently, no extra options are defined.
    @rtype: L{ParsedDocstring}
    """
    writer = _DocumentPseudoWriter()
    reader = _EpydocReader(errors) # Outputs errors to the list.
    publish_string(docstring, writer=writer, reader=reader)
    return ParsedRstDocstring(writer.document)

class ParsedRstDocstring(ParsedDocstring):
    """
    An encoded version of a ReStructuredText docstring.  The contents
    of the docstring are encoded in the L{_document} instance
    variable.

    @ivar _document: A ReStructuredText document, encoding the
        docstring.
    @type _document: L{docutils.nodes.document}
    """
    def __init__(self, document):
        """
        @type document: L{docutils.nodes.document}
        """
        self._document = document

    def split_fields(self, errors=None):
        # Inherit docs
        visitor = _SplitFieldsTranslator(self._document, errors)
        self._document.walk(visitor)
        return self, visitor.fields

    def summary(self):
        # Inherit docs
        visitor = _SummaryExtractor(self._document)
        self._document.walk(visitor)
        return visitor.summary

#     def concatenate(self, other):
#         result = self._document.copy()
#         for child in self._document.children + other._document.children:
#             visitor = TreeCopyVisitor(self._document)
#             child.walkabout(visitor)
#             result.append(visitor.get_tree_copy())
#         return ParsedRstDocstring(result)
        
    def to_html(self, docstring_linker, **options):
        # Inherit docs
        visitor = _EpydocHTMLTranslator(self._document, docstring_linker)
        self._document.walkabout(visitor)
        return ''.join(visitor.body)

    def to_latex(self, docstring_linker, **options):
        # Inherit docs
        visitor = _EpydocLaTeXTranslator(self._document, docstring_linker)
        self._document.walkabout(visitor)
        return ''.join(visitor.body)

    def to_plaintext(self, docstring_linker, **options):
        # This is should be replaced by something better:
        return self._document.astext() 

    def __repr__(self): return '<ParsedRstDocstring: ...>'

class _EpydocReader(StandaloneReader):
    """
    A reader that captures all errors that are generated by parsing,
    and appends them to a list.
    """
    def __init__(self, errors):
        self._errors = errors
        StandaloneReader.__init__(self)
        
    def new_document(self):
        document = new_document(self.source.source_path, self.settings)
        document.reporter.attach_observer(self.report)
        document.reporter.set_conditions('', 10000, 10000, None)
        self._encoding = document.reporter.encoding
        self._error_handler = document.reporter.error_handler
        return document

    def report(self, error):
        try: is_fatal = int(error['level']) > 2
        except: is_fatal = 1
        try: linenum = int(error['line'])
        except: linenum = None

        msg = ''.join([c.astext().encode(self._encoding, self._error_handler)
                       for c in error.children])

        self._errors.append(ParseError(msg, linenum, is_fatal))
        
class _DocumentPseudoWriter(Writer):
    """
    A pseudo-writer for the docutils framework, that can be used to
    access the document itself.  The output of C{_DocumentPseudoWriter}
    is just an empty string; but after it has been used, the most
    recently processed document is available as the instance variable
    C{document}

    @type document: L{docutils.nodes.document}
    @ivar document: The most recently processed document.
    """
    def __init__(self):
        self.document = None
        Writer.__init__(self)
        
    def translate(self):
        self.output = ''
        
class _SummaryExtractor(NodeVisitor):
    """
    A docutils node visitor that extracts the first sentence from
    the first paragraph in a document.
    """
    def __init__(self, document):
        NodeVisitor.__init__(self, document)
        self.summary = None
        
    def visit_document(self, node):
        self.summary = None
        
    def visit_paragraph(self, node):
        if self.summary is not None: return

        summary_pieces = []
        # Extract the first sentence.
        for child in node.children:
            if isinstance(child, docutils.nodes.Text):
                m = re.match(r'(\s*[\w\W]*?\.)(\s|$)', child.data)
                if m:
                    summary_pieces.append(docutils.nodes.Text(m.group(1)))
                    break
            summary_pieces.append(child)
            
        summary_doc = self.document.copy()
        summary_doc.children = summary_pieces
        self.summary = ParsedRstDocstring(summary_doc)

    def unknown_visit(self, node):
        'Ignore all unknown nodes'

class _SplitFieldsTranslator(NodeVisitor):
    """
    A docutils translator that removes all fields from a document, and
    collects them into the instance variable C{fields}

    @ivar fields: The fields of the most recently walked document.
    @type fields: C{list} of L{Field<markup.Field>}
    """
    def __init__(self, document, errors):
        NodeVisitor.__init__(self, document)
        self._errors = errors
        self.fields = []
        self._newfields = {}

    def visit_document(self, node):
        self.fields = []

    def visit_field(self, node):
        # Remove the field from the tree.
        node.parent.remove(node)

        # Extract the field name & optional argument
        tag = node.children[0].astext().split(None, 1)
        tagname = tag[0]
        if len(tag)>1: arg = tag[1]
        else: arg = None

        # Handle special fields:
        fbody = node.children[1].children
        if arg is None:
            for (list_tag, entry_tag) in CONSOLIDATED_FIELDS.items():
                if tagname.lower() == list_tag:
                    try:
                        self.handle_consolidated_field(fbody, entry_tag)
                        return
                    except ValueError, e:
                        estr = 'Unable to split consolidated field '
                        estr += '"%s" - %s' % (tagname, e)
                        self._errors.append(ParseError(estr, is_fatal=0))
                        
                        # Use a @newfield to let it be displayed as-is.
                        if not self._newfields.has_key(tagname.lower()):
                            newfield = Field('newfield', tagname.lower(),
                                             parse(tagname, 'plaintext'))
                            self.fields.append(newfield)
                            self._newfields[tagname.lower()] = 1
                        
        self._add_field(tagname, arg, fbody)

    def _add_field(self, tagname, arg, fbody):
        field_doc = self.document.copy()
        for child in fbody: field_doc.append(child)
        field_pdoc = ParsedRstDocstring(field_doc)
        self.fields.append(Field(tagname, arg, field_pdoc))
            
    def visit_field_list(self, node):
        # Remove the field list from the tree.  The visitor will still walk
        # over the node's children.
        node.parent.remove(node)

    def handle_consolidated_field(self, body, tagname):
        """
        Attempt to handle a consolidated section.
        """
        # Check that it contains a bulleted list.
        if len(body) != 1 or body[0].tagname != 'bullet_list':
            raise ValueError('field does not contain a single bulleted list.')

        # Check that each list item begins with interpreted text
        n = 0
        for item in body[0].children:
            n += 1
            if item.tagname != 'list_item':
                raise ValueError('bad bulleted list (bad child).')
            if len(item.children) == 0: 
                raise ValueError('bad bulleted list (empty).')
            if item.children[0].tagname != 'paragraph':
                if item.children[0].tagname == 'definition_list':
                    raise ValueError(('list item %d contains a definition '+
                                      'list (it\'s probably indented '+
                                      'wrong).') % n)
                else:
                    raise ValueError(('list item %d does not begin with '+
                                      'an identifier.') % n)
            if len(item.children[0].children) == 0: 
                raise ValueError(('list item %d does not begin with '+
                                  'an identifier.') % n)
            if item.children[0].children[0].tagname != 'title_reference':
                raise ValueError(('list item %d does not begin with '+
                                  'an identifier.') % n)

        # Everything looks good; convert to multiple fields.
        for item in body[0].children:
            # Extract the arg
            arg = item.children[0].children[0].astext()

            # Extract the field body, and remove the arg
            fbody = item.children[:]
            fbody[0] = fbody[0].copy()
            fbody[0].children = item.children[0].children[1:]

            # Remove the separating ":", if present
            if (len(fbody[0].children) > 0 and
                isinstance(fbody[0].children[0], docutils.nodes.Text)):
                child = fbody[0].children[0]
                if child.data[:1] in ':-':
                    child.data = child.data[1:].lstrip()
                elif child.data[:2] == ' -':
                    child.data = child.data[2:].lstrip()

            # Wrap the field body, and add a new field
            self._add_field(tagname, arg, fbody)
        
    def unknown_visit(self, node):
        'Ignore all unknown nodes'

class _EpydocLaTeXTranslator(LaTeXTranslator):
    def __init__(self, document, docstring_linker):
        # Set the document's settings.
        settings = OptionParser([LaTeXWriter()]).get_default_values()
        document.settings = settings

        LaTeXTranslator.__init__(self, document)
        self._linker = docstring_linker

        # Start at section level 3.
        self.section_level = 3

    # Handle interpreted text (crossreferences)
    def visit_title_reference(self, node):
        target = self.encode(node.astext())
        xref = self._linker.translate_identifier_xref(target, target)
        self.body.append(xref)
        raise SkipNode

    def visit_document(self, node): pass
    def depart_document(self, node): pass
        
class _EpydocHTMLTranslator(HTMLTranslator):
    def __init__(self, document, docstring_linker):
        # Set the document's settings.
        settings = OptionParser([HTMLWriter()]).get_default_values()
        document.settings = settings
    
        HTMLTranslator.__init__(self, document)
        self._linker = docstring_linker

    # Handle interpreted text (crossreferences)
    def visit_title_reference(self, node):
        target = self.encode(node.astext())
        xref = self._linker.translate_identifier_xref(target, target)
        self.body.append(xref)
        raise SkipNode

    def visit_document(self, node): pass
    def depart_document(self, node): pass
        
    def starttag(self, node, tagname, suffix='\n', infix='', **attributes):
        """
        This modified version of starttag makes a few changes to HTML
        tags, to prevent them from conflicting with epydoc.  In particular:
          - existing class attributes are prefixed with C{'rst-'}
          - existing names are prefixed with C{'rst-'}
          - hrefs starting with C{'#'} are prefixed with C{'rst-'}
          - all headings (C{<hM{n}>}) are given the css class C{'heading'}
        """
        # Prefix all CSS classes with "rst-"
        if attributes.has_key('class'):
            attributes['class'] = 'rst-%s' % attributes['class']

        # Prefix all names with "rst-", to avoid conflicts
        if attributes.has_key('id'):
            attributes['id'] = 'rst-%s' % attributes['id']
        if attributes.has_key('name'):
            attributes['name'] = 'rst-%s' % attributes['name']
        if attributes.has_key('href') and attributes['href'][:1]=='#':
            attributes['href'] = '#rst-%s' % attributes['href'][1:]

        # For headings, use class="heading"
        if re.match(r'^h\d+$', tagname):
            attributes['class'] = 'heading'
        
        return HTMLTranslator.starttag(self, node, tagname, suffix,
                                       infix, **attributes)
