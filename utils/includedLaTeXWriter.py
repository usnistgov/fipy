## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "includedLaTeXWriter.py"
 #                                    created: 9/29/04 {11:38:07 AM} 
 #                                last update: 10/27/04 {2:45:00 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James Warren <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  includedLaTeXWriter.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-09-29 JEG 1.0 original
 # ###################################################################
 ##

from docutils.writers.latex2e import LaTeXTranslator, Writer as LaTeXWriter
from docutils import languages

class NotStupidLaTeXTranslator(LaTeXTranslator):
    import re
    externalNISThref = re.compile(r"http://www\.nist\.gov/cgi-bin/exit_nist\.cgi\?url=(.*)")
    
    def visit_reference(self, node):
	# BUG: hash_char "#" is trouble some in LaTeX.
	# mbox and other environment do not like the '#'.
	hash_char = '\\#'
	if node.has_key('refuri'):
	    href = node['refuri'].replace('#',hash_char)
	elif node.has_key('refid'):
	    href = hash_char + node['refid']
	elif node.has_key('refname'):
	    href = hash_char + self.document.nameids[node['refname']]
	else:
	    raise AssertionError('Unknown reference.')
	
	# It doesn't make any sense to redirect external 
	# web pages through a NIST disclaimer page
	# when linked from a PDF document.
	externalMatch = self.externalNISThref.match(href)
	if externalMatch:
	    href = externalMatch.group(1)
	    
	self.body.append('\\href{%s}{' % href)

class IncludedLaTeXWriter(LaTeXWriter):
    def write(self, document, destination):
	self.document = document
	self.language = languages.get_language(
	    document.settings.language_code)
	self.destination = destination
	self.translate()
	output = self.destination.write(''.join(self.body))
	return output
	
    def translate(self):
	visitor = NotStupidLaTeXTranslator(self.document)
	self.document.walkabout(visitor)
	self.output = visitor.astext()
	self.head_prefix = visitor.head_prefix
	self.head = visitor.head
	self.body_prefix = visitor.body_prefix
	self.body = visitor.body
	self.body_suffix = visitor.body_suffix

