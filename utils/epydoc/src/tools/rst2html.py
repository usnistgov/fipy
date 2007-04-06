#!/usr/bin/python

r"""
A customized driver for converting docutils reStructuredText documents
into HTML.  This is used to generated HTML versions of the regression
files, for the webpage.
"""

# Docutils imports
from docutils.core import publish_cmdline, default_description
from docutils.writers.html4css1 import HTMLTranslator, Writer as HTMLWriter
import docutils.nodes

# Epydoc imports.  Make sure path contains the 'right' epydoc.
import sys, os
sys.path.insert(0, '../')
import epydoc.markup.restructuredtext   # register the 'python' directive
from epydoc.markup.doctest import doctest_to_html, doctest_to_latex, \
                                  HTMLDoctestColorizer

from epydoc.docwriter.xlink import ApiLinkReader

class CustomizedReader(ApiLinkReader):
    settings_defaults = (ApiLinkReader.settings_defaults or {}).copy()
    settings_defaults.update({
        'external_api': [ 'epydoc' ],
        'external_api_root': [ 'epydoc:http://epydoc.sourceforge.net/api/' ],
        'external_api_file': [ 'epydoc:' + os.path.join(
            os.path.split(__file__)[0], '../../html/api/api-objects.txt') ],
        })
    
class CustomizedHTMLWriter(HTMLWriter):
    settings_defaults = (HTMLWriter.settings_defaults or {}).copy()
    settings_defaults.update({
        'stylesheet': 'custom.css',
        'stylesheet_path': None,
        'output_encoding': 'ascii',
        'output_encoding_error_handler': 'xmlcharrefreplace',
        'embed_stylesheet': False,
        })
        
    def __init__(self):
        HTMLWriter.__init__(self)
        self.translator_class = CustomizedHTMLTranslator

class CustomizedHTMLTranslator(HTMLTranslator):
    def visit_doctest_block(self, node):
        pysrc = node[0].astext()
        if node.get('codeblock'):
            self.body.append(HTMLDoctestColorizer().colorize_codeblock(pysrc))
        else:
            self.body.append(doctest_to_html(pysrc))
        raise docutils.nodes.SkipNode

description = ('Generates HTML documents from reStructuredText '
               'documents.  ' + default_description)
writer = CustomizedHTMLWriter()
reader = CustomizedReader()

#this doesn't work. Put ``.. default-role:: epydoc`` in the doctests instead.
#docutils.parsers.rst.roles.DEFAULT_INTERPRETED_ROLE = 'epydoc'

docutils.core.publish_cmdline(reader=reader, writer=writer,
                              description=description)
