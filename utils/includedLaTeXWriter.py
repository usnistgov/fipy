## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "includedLaTeXWriter.py"
 #                                    created: 9/29/04 {11:38:07 AM} 
 #                                last update: 9/15/05 {6:42:58 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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
from docutils import languages, nodes

class NotStupidLaTeXTranslator(LaTeXTranslator):
    import re
    externalNISThref = re.compile(r"http://www\.nist\.gov/cgi-bin/exit_nist\.cgi\?url=(.*)")
    
    def __init__(self, document):
        LaTeXTranslator.__init__(self, document)
        self.d_class._class_sections['startlower'] = ('subsection', 'subsubsection')
        
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
        
    def visit_image(self, node):
        """
        reST isn't smart enough to direct different output formats to use different 
        image file extensions: http://thread.gmane.org/gmane.text.docutils.user/1239
        
        For LaTeX output, we just strip off the file appendage, as the graphicx 
        package is smart enough.
        """
        attrs = node.attributes
        # Add image URI to dependency list, assuming that it's
        # referring to a local file.
        self.settings.record_dependencies.add(attrs['uri'])
        import os
        imageURI = os.path.splitext(attrs['uri'])[0]
        print attrs['uri'], imageURI
        pre = []                        # in reverse order
        post = ['\\includegraphics{%s}' % imageURI]
        inline = isinstance(node.parent, nodes.TextElement)
        if attrs.has_key('scale'):
            # Could also be done with ``scale`` option to
            # ``\includegraphics``; doing it this way for consistency.
            pre.append('\\scalebox{%f}{' % (attrs['scale'] / 100.0,))
            post.append('}')
        if attrs.has_key('align'):
            align_prepost = {
                # By default latex aligns the top of an image.
                (1, 'top'): ('', ''),
                (1, 'middle'): ('\\raisebox{-0.5\\height}{', '}'),
                (1, 'bottom'): ('\\raisebox{-\\height}{', '}'),
                (0, 'center'): ('{\\hfill', '\\hfill}'),
                # These 2 don't exactly do the right thing.  The image should
                # be floated alongside the paragraph.  See
                # http://www.w3.org/TR/html4/struct/objects.html#adef-align-IMG
                (0, 'left'): ('{', '\\hfill}'),
                (0, 'right'): ('{\\hfill', '}'),}
            try:
                pre.append(align_prepost[inline, attrs['align']][0])
                post.append(align_prepost[inline, attrs['align']][1])
            except KeyError:
                pass                    # XXX complain here?
        if not inline:
            pre.append('\n')
            post.append('\n')
        pre.reverse()
        self.body.extend(pre + post)

    def visit_admonition(self, node, name=''):
        self.body.append('\\begin{reSTadmonition}\n')
        if name:
            self.body.append('['+ self.language.labels[name] + ']');
        self.body.append('\n')


    def depart_admonition(self, node=None):
        self.body.append('\\end{reSTadmonition}\n');


class IncludedLaTeXWriter(LaTeXWriter):
    def __init__(self):
        LaTeXWriter.__init__(self)
        self.translator_class = NotStupidLaTeXTranslator

    def translate(self):
        LaTeXWriter.translate(self)
        self.output = ''.join(self.body)


