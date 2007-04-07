#
# epydoc.py: epydoc LaTeX output generator
# Edward Loper
#
# Created [01/30/01 05:18 PM]
# $Id: latex.py 1587 2007-03-20 00:58:51Z dvarrazzo $
#

"""
The LaTeX output generator for epydoc.  The main interface provided by
this module is the L{LatexWriter} class.

@todo: Inheritance=listed
"""
__docformat__ = 'epytext en'

import os.path, sys, time, re, textwrap, codecs

from epydoc.apidoc import *
from epydoc.compat import *
import epydoc
from epydoc import log
from epydoc import markup
from epydoc.util import plaintext_to_latex
import epydoc.markup

class LatexWriter:
    PREAMBLE = [
        "\\documentclass{article}",
        "\\usepackage[%s]{epydoc}",
        ]

    SECTIONS = ['\\part{%s}', '\\chapter{%s}', '\\section{%s}',
                '\\subsection{%s}', '\\subsubsection{%s}',
                '\\textbf{%s}']

    STAR_SECTIONS = ['\\part*{%s}', '\\chapter*{%s}', '\\section*{%s}',
                     '\\subsection*{%s}', '\\subsubsection*{%s}',
                     '\\textbf{%s}']

    def __init__(self, docindex, **kwargs):
        self.docindex = docindex
        # Process keyword arguments
        self._show_private = kwargs.get('private', 0)
        self._prj_name = kwargs.get('prj_name', None) or 'API Documentation'
        self._crossref = kwargs.get('crossref', 1)
        self._index = kwargs.get('index', 1)
        self._list_classes_separately=kwargs.get('list_classes_separately',0)
        self._inheritance = kwargs.get('inheritance', 'listed')
        self._exclude = kwargs.get('exclude', 1)
        self._list_submodules = kwargs.get('list_submodules', 1)
        self._top_section = 2
        self._index_functions = 1
        self._hyperref = 1
        self._pdflatex = kwargs.get('pdflatex', 0)

        #: The Python representation of the encoding.
        #: Update L{latex_encodings} in case of mismatch between it and
        #: the C{inputenc} LaTeX package.
        self._encoding = kwargs.get('encoding', 'utf-8')

        self.valdocs = sorted(docindex.reachable_valdocs(
            imports=False, packages=False, bases=False, submodules=False, 
            subclasses=False, private=self._show_private))
        self._num_files = self.num_files()
        # For use with select_variables():
        if self._show_private: self._public_filter = None
        else: self._public_filter = True

    def write(self, directory=None):
        """
        Write the API documentation for the entire project to the
        given directory.

        @type directory: C{string}
        @param directory: The directory to which output should be
            written.  If no directory is specified, output will be
            written to the current directory.  If the directory does
            not exist, it will be created.
        @rtype: C{None}
        @raise OSError: If C{directory} cannot be created,
        @raise OSError: If any file cannot be created or written to.
        """
        # For progress reporting:
        self._files_written = 0.
        
        # Set the default values for ValueDoc formatted representations.
        orig_valdoc_defaults = (ValueDoc.SUMMARY_REPR_LINELEN,
                                ValueDoc.REPR_LINELEN,
                                ValueDoc.REPR_MAXLINES)
        ValueDoc.SUMMARY_REPR_LINELEN = 60
        ValueDoc.REPR_LINELEN = 52
        ValueDoc.REPR_MAXLINES = 5

        # Create destination directories, if necessary
        if not directory: directory = os.curdir
        self._mkdir(directory)
        self._directory = directory
        
        # Write the top-level file.
        self._write(self.write_topfile, directory, 'api.tex')

        # Write the module & class files.
        for val_doc in self.valdocs:
            if isinstance(val_doc, ModuleDoc):
                filename = '%s-module.tex' % val_doc.canonical_name
                self._write(self.write_module, directory, filename, val_doc)
            elif (isinstance(val_doc, ClassDoc) and 
                  self._list_classes_separately):
                filename = '%s-class.tex' % val_doc.canonical_name
                self._write(self.write_class, directory, filename, val_doc)

        # Restore defaults that we changed.
        (ValueDoc.SUMMARY_REPR_LINELEN, ValueDoc.REPR_LINELEN,
         ValueDoc.REPR_MAXLINES) = orig_valdoc_defaults
        
    def _write(self, write_func, directory, filename, *args):
        # Display our progress.
        self._files_written += 1
        log.progress(self._files_written/self._num_files, filename)
        
        path = os.path.join(directory, filename)
        if self._encoding == 'utf-8':
            f = codecs.open(path, 'w', 'utf-8')
            write_func(f.write, *args)
            f.close()
        else:
            result = []
            write_func(result.append, *args)
            s = u''.join(result)
            try:
                s = s.encode(self._encoding)
            except UnicodeError:
                log.error("Output could not be represented with the "
                          "given encoding (%r).  Unencodable characters "
                          "will be displayed as '?'.  It is recommended "
                          "that you use a different output encoding (utf-8, "
                          "if it's supported by latex on your system)."
                          % self._encoding)
                s = s.encode(self._encoding, 'replace')
            f = open(path, 'w')
            f.write(s)
            f.close()

    def num_files(self):
        """
        @return: The number of files that this C{LatexFormatter} will
            generate.
        @rtype: C{int}
        """
        n = 1
        for doc in self.valdocs:
            if isinstance(doc, ModuleDoc): n += 1
            if isinstance(doc, ClassDoc) and self._list_classes_separately:
                n += 1
        return n
        
    def _mkdir(self, directory):
        """
        If the given directory does not exist, then attempt to create it.
        @rtype: C{None}
        """
        if not os.path.isdir(directory):
            if os.path.exists(directory):
                raise OSError('%r is not a directory' % directory)
            os.mkdir(directory)
            
    #////////////////////////////////////////////////////////////
    #{ Main Doc File
    #////////////////////////////////////////////////////////////

    def write_topfile(self, out):
        self.write_header(out, 'Include File')
        self.write_preamble(out)
        out('\n\\begin{document}\n\n')
        self.write_start_of(out, 'Header')

        # Write the title.
        self.write_start_of(out, 'Title')
        out('\\title{%s}\n' % plaintext_to_latex(self._prj_name, 1))
        out('\\author{API Documentation}\n')
        out('\\maketitle\n')

        # Add a table of contents.
        self.write_start_of(out, 'Table of Contents')
        out('\\addtolength{\\parskip}{-1ex}\n')
        out('\\tableofcontents\n')
        out('\\addtolength{\\parskip}{1ex}\n')

        # Include documentation files.
        self.write_start_of(out, 'Includes')
        for val_doc in self.valdocs:
            if isinstance(val_doc, ModuleDoc):
                out('\\include{%s-module}\n' % val_doc.canonical_name)

        # If we're listing classes separately, put them after all the
        # modules.
        if self._list_classes_separately:
            for val_doc in self.valdocs:
                if isinstance(val_doc, ClassDoc):
                    out('\\include{%s-class}\n' % val_doc.canonical_name)

        # Add the index, if requested.
        if self._index:
            self.write_start_of(out, 'Index')
            out('\\printindex\n\n')

        # Add the footer.
        self.write_start_of(out, 'Footer')
        out('\\end{document}\n\n')

    def write_preamble(self, out):
        # If we're generating an index, add it to the preamble.
        if self._index:
            out('\n'.join(self.PREAMBLE) % 'index')
        else:
            out('\n'.join(self.PREAMBLE) % '')    
        out('\n')
        
        # Set the encoding.
        out('\\usepackage[%s]{inputenc}' % self.get_latex_encoding())

        # If we're generating hyperrefs, add the appropriate packages.
        # !!!!!!!!!!!!!!!!!!!!!!
        # !!! JEG - this needs to be the last thing in the preamble
        # !!!!!!!!!!!!!!!!!!!!!!
        if self._hyperref:
            out('\\definecolor{UrlColor}{rgb}{0,0.08,0.45}\n')
            
            if self._pdflatex:
                driver = 'pdftex'
            else:
                driver = 'dvips'
                
            out('\\usepackage[%s, pagebackref, pdftitle={%s}, '
                'pdfcreator={epydoc %s}, bookmarks=true, '
                'bookmarksopen=false, pdfpagemode=UseOutlines, '
                'colorlinks=true, linkcolor=black, anchorcolor=black, '
                'citecolor=black, filecolor=black, menucolor=black, '
                'pagecolor=black, urlcolor=UrlColor]{hyperref}\n' %
                (driver, self._prj_name or '', epydoc.__version__))
            
        # If restructuredtext was used, then we need to extend
        # the prefix to include LatexTranslator.head_prefix.
        if 'restructuredtext' in epydoc.markup.MARKUP_LANGUAGES_USED:
            from epydoc.markup import restructuredtext
            rst_head = restructuredtext.latex_head_prefix()
            for line in rst_head[1:]:
                m = re.match(r'\\usepackage(\[.*?\])?{(.*?)}', line)
                if m and m.group(2) in (
                    'babel', 'hyperref', 'color', 'alltt', 'parskip',
                    'fancyhdr', 'boxedminipage', 'makeidx',
                    'multirow', 'longtable', 'tocbind', 'assymb',
                    'fullpage', 'inputenc'):
                    pass
                else:
                    out(line)

        
    #////////////////////////////////////////////////////////////
    #{ Chapters
    #////////////////////////////////////////////////////////////

    def write_module(self, out, doc):
        self.write_header(out, doc)
        self.write_start_of(out, 'Module Description')

        # Add this module to the index.
        out('    ' + self.indexterm(doc, 'start'))

        # Add a section marker.
        out(self.section('%s %s' % (self.doc_kind(doc),
                                    self._dotted(doc.canonical_name)),
                         ref=doc))

        # Add the module's description.
        if doc.descr not in (None, UNKNOWN):
            out(self.docstring_to_latex(doc.descr))

        # Add version, author, warnings, requirements, notes, etc.
        self.write_standard_fields(out, doc)

        # If it's a package, list the sub-modules.
        if self._list_submodules and doc.submodules != UNKNOWN and doc.submodules:
            self.write_module_list(out, doc)

        # Contents.
        if self._list_classes_separately:
            self.write_class_list(out, doc)
        self.write_func_list(out, 'Functions', doc, 'function')
        self.write_var_list(out, 'Variables', doc, 'other')

        # Class list.
        if not self._list_classes_separately:
            classes = doc.select_variables(imported=False, value_type='class',
                                           public=self._public_filter)
            for var_doc in classes:
                self.write_class(out, var_doc.value)

        # Mark the end of the module (for the index)
        out('    ' + self.indexterm(doc, 'end'))

    def write_class(self, out, doc):
        if self._list_classes_separately:
            self.write_header(out, doc)
        self.write_start_of(out, 'Class Description')

        # Add this class to the index.
        out('    ' + self.indexterm(doc, 'start'))

        # Add a section marker.
        if self._list_classes_separately:
            seclevel = 0
            out(self.section('%s %s' % (self.doc_kind(doc),
                                        self._dotted(doc.canonical_name)), 
                             seclevel, ref=doc))
        else:
            seclevel = 1
            out(self.section('%s %s' % (self.doc_kind(doc),
                                        self._dotted(doc.canonical_name[-1])), 
                             seclevel, ref=doc))

        # Add our base list.
        if doc.bases not in (UNKNOWN, None) and len(doc.bases) > 0:
            out(self.base_tree(doc))

        # The class's known subclasses
        if doc.subclasses not in (UNKNOWN, None) and len(doc.subclasses) > 0:
            sc_items = [self._hyperlink(sc, '%s' % sc.canonical_name)
                        for sc in doc.subclasses]
            out(self._descrlist(sc_items, 'Known Subclasses', short=1))

        # The class's description.
        if doc.descr not in (None, UNKNOWN):
            out(self.docstring_to_latex(doc.descr))

        # Version, author, warnings, requirements, notes, etc.
        self.write_standard_fields(out, doc)

        # Contents.
        self.write_func_list(out, 'Methods', doc, 'method',
                             seclevel+1)
        self.write_var_list(out, 'Properties', doc,
                            'property', seclevel+1)
        self.write_var_list(out, 'Class Variables', doc, 
                            'classvariable', seclevel+1)
        self.write_var_list(out, 'Instance Variables', doc, 
                            'instancevariable', seclevel+1)

        # Mark the end of the class (for the index)
        out('    ' + self.indexterm(doc, 'end'))

    #////////////////////////////////////////////////////////////
    #{ Module hierarchy trees
    #////////////////////////////////////////////////////////////
    
    def write_module_tree(self, out):
        modules = [doc for doc in self.valdocs
                   if isinstance(doc, ModuleDoc)]
        if not modules: return
        
        # Write entries for all top-level modules/packages.
        out('\\begin{itemize}\n')
        out('\\setlength{\\parskip}{0ex}\n')
        for doc in modules:
            if (doc.package in (None, UNKNOWN) or
                doc.package not in self.valdocs):
                self.write_module_tree_item(out, doc)
        return s +'\\end{itemize}\n'

    def write_module_list(self, out, doc):
        if len(doc.submodules) == 0: return
        self.write_start_of(out, 'Modules')
        
        out(self.section('Modules', 1))
        out('\\begin{EpydocModuleList}\n')

        for group_name in doc.group_names():
            if not doc.submodule_groups[group_name]: continue
            if group_name:
                out('  \\item %s\n' % group_name)
                out('  \\begin{EpydocModuleSubList}\n')
            for submodule in doc.submodule_groups[group_name]:
                self.write_module_tree_item(out, submodule)
            if group_name:
                out('  \end{EpydocModuleSubList}\n')

        out('\\end{EpydocModuleList}\n\n')

    def write_module_tree_item(self, out, doc, depth=0):
        """
        Helper function for L{write_module_tree} and L{write_module_list}.
        
        @rtype: C{string}
        """
        out(' '*depth + '\\item[%s]' % self._hyperlink(doc, doc.canonical_name[-1]))

        if doc.summary not in (None, UNKNOWN):
            out(' %s\n' % self.docstring_to_latex(doc.summary))
        if self._crossref:
            out('\\CrossRef{%s}\n\n' % self.label(doc))
        if doc.submodules != UNKNOWN and doc.submodules:
            out(' '*depth + '  \\begin{EpydocModuleSubList}\n')
            for submodule in doc.submodules:
                self.write_module_tree_item(out, submodule, depth+4)
            out(' '*depth + '  \\end{EpydocModuleSubList}\n')

    #////////////////////////////////////////////////////////////
    #{ Base class trees
    #////////////////////////////////////////////////////////////

    def base_tree(self, doc, width=None, linespec=None):
        if width is None:
            width = self._find_tree_width(doc)+2
            linespec = []
            s = ('&'*(width-4)+'\\multicolumn{2}{l}{\\textbf{%s}}\n' %
                   self._dotted('%s'%self._base_name(doc)))
            s += '\\end{tabular}\n\n'
            top = 1
        else:
            s = self._base_tree_line(doc, width, linespec)
            top = 0
        
        if isinstance(doc, ClassDoc):
            for i in range(len(doc.bases)-1, -1, -1):
                base = doc.bases[i]
                spec = (i > 0)
                s = self.base_tree(base, width, [spec]+linespec) + s

        if top:
            s = '\\begin{tabular}{%s}\n' % (width*'c') + s

        return s

    def _base_name(self, doc):
        if doc.canonical_name is None:
            if doc.parse_repr is not None:
                return doc.parse_repr
            else:
                return '??'
        else:
            return '%s' % doc.canonical_name

    def _find_tree_width(self, doc):
        if not isinstance(doc, ClassDoc): return 2
        width = 2
        for base in doc.bases:
            width = max(width, self._find_tree_width(base)+2)
        return width

    def _base_tree_line(self, doc, width, linespec):
        base_name = self._dotted(self._base_name(doc))
        
        # linespec is a list of booleans.
        s = '%% Line for %s, linespec=%s\n' % (base_name, linespec)

        labelwidth = width-2*len(linespec)-2

        # The base class name.
        s += ('\\multicolumn{%s}{r}{' % labelwidth)
        s += '\\settowidth{\\BCL}{%s}' % base_name
        s += '\\multirow{2}{\\BCL}{%s}}\n' % self._hyperlink(doc, self._base_name(doc))

        # The vertical bars for other base classes (top half)
        for vbar in linespec:
            if vbar: s += '&&\\multicolumn{1}{|c}{}\n'
            else: s += '&&\n'

        # The horizontal line.
        s += '  \\\\\\cline{%s-%s}\n' % (labelwidth+1, labelwidth+1)

        # The vertical bar for this base class.
        s += '  ' + '&'*labelwidth
        s += '\\multicolumn{1}{c|}{}\n'

        # The vertical bars for other base classes (bottom half)
        for vbar in linespec:
            if vbar: s += '&\\multicolumn{1}{|c}{}&\n'
            else: s += '&&\n'
        s += '  \\\\\n'

        return s
        
    #////////////////////////////////////////////////////////////
    #{ Class List
    #////////////////////////////////////////////////////////////
    
    def write_class_list(self, out, doc):
        groups = [(plaintext_to_latex(group_name),
                   doc.select_variables(group=group_name, imported=False,
                                        value_type='class',
                                        public=self._public_filter))
                  for group_name in doc.group_names()]

        # Discard any empty groups; and return if they're all empty.
        groups = [(g,vars) for (g,vars) in groups if vars]
        if not groups: return

        # Write a header.
        self.write_start_of(out, 'Classes')
        out(self.section('Classes', 1))
        out('\\begin{EpydocClassList}\n')

        for name, var_docs in groups:
            if name:
                out('  \\item[%s]\n' % name)
                out('  \\begin{EpydocClassSubList}\n')
            # Add the lines for each class
            for var_doc in var_docs:
                self.write_class_list_line(out, var_doc)
            if name:
                out('  \\end{EpydocClassSubList}\n')

        out('\\end{EpydocClassList}\n')

    def write_class_list_line(self, out, var_doc):
        if var_doc.value in (None, UNKNOWN): return # shouldn't happen
        doc = var_doc.value
        out('  ' + '\\item[%s]' % self._hyperlink(var_doc.target, 
                                                  var_doc.name))
        if doc.summary not in (None, UNKNOWN):
            out(': %s\n' % self.docstring_to_latex(doc.summary))
        if self._crossref:
            out(('\CrossRef{%s}\n\n' % self.label(doc)))
        
    #////////////////////////////////////////////////////////////
    #{ Function List
    #////////////////////////////////////////////////////////////
    
    def write_func_list(self, out, heading, doc, value_type, seclevel=1):
        groups = [(plaintext_to_latex(group_name),
                   doc.select_variables(group=group_name, imported=False,
                                        value_type=value_type,
                                        public=self._public_filter))
                  for group_name in doc.group_names()]

        # Discard any empty groups; and return if they're all empty.
        groups = [(g,vars) for (g,vars) in groups if vars]
        if not groups: return

        # Write a header.
        self.write_start_of(out, heading)
        out('  '+self.section(heading, seclevel))

        for name, var_docs in groups:
            if name:
                out('\\begin{EpydocFunctionGroup}[%s]\n' % (name or ''))
            for var_doc in var_docs:
                self.write_func_list_box(out, var_doc)
            # [xx] deal with inherited methods better????
            #if (self._inheritance == 'listed' and
            #    isinstance(container, ClassDoc)):
            #    out(self._inheritance_list(group, container.uid()))
            if name:
                out('\\end{EpydocFunctionGroup}\n')

            
    def write_func_list_box(self, out, var_doc):
        func_doc = var_doc.value
        is_inherited = (var_doc.overrides not in (None, UNKNOWN))

        out('\\begin{EpydocFunction}')
        # Function signature.
        out('%s\n\n' % self.function_signature(var_doc))

        # nb: this gives the containing section, not a reference
        # directly to the function.
        if not is_inherited:
            out('    %s\n' % self.indexterm(func_doc))

        # Description
        if func_doc.descr not in (None, UNKNOWN):
            out(' '*4 + '\\begin{EpydocDocstring}\n')
            out(self.docstring_to_latex(func_doc.descr, 4))
            out(' '*4 + '\\end{EpydocDocstring}\n')

        # Parameters
        if func_doc.arg_descrs or func_doc.arg_types:
            # Find the longest name.
            longest = max([0]+[len(n) for n in func_doc.arg_types])
            for names, descrs in func_doc.arg_descrs:
                longest = max([longest]+[len(n) for n in names])
            # Table header.
            out(' '*6+'\\begin{EpydocFunctionParameters}{%s}\n' % (longest*'x'))
            # Params that have @type but not @param info:
            unseen_types = set(func_doc.arg_types)
            # List everything that has a @param:
            for (arg_names, arg_descr) in func_doc.arg_descrs:
                arg_name = plaintext_to_latex(', '.join(arg_names))
                out('%s\\item[%s]\n\n' % (' '*10, arg_name))
                out(self.docstring_to_latex(arg_descr, 10))
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # !!! JEG - this loop needs abstracting
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                for arg_name in arg_names:
                    arg_typ = func_doc.arg_types.get(arg_name)
                    if arg_typ is not None:
                        if len(arg_names) == 1:
                            lhs = 'type'
                        else:
                            lhs = 'type of %s' % arg_name
                        rhs = self.docstring_to_latex(arg_typ).strip()
                        out('%s\\textit{(%s=%s)}\n\n' % (' '*12, lhs, rhs))
            out(' '*6+'\\end{EpydocFunctionParameters}\n\n')
                
        # Returns
        rdescr = func_doc.return_descr
        rtype = func_doc.return_type
        if rdescr not in (None, UNKNOWN) or rtype not in (None, UNKNOWN):
            out(' '*6+'\\EpydocFunctionReturns')
            if rtype not in (None, UNKNOWN):
                out('[%s]' % self.docstring_to_latex(rtype, 6).strip())
            if rdescr not in (None, UNKNOWN):
                out('{%s}' % self.docstring_to_latex(rdescr, 6))
            else:
                out('{}')
            out('\n\n')

        # Raises
        if func_doc.exception_descrs not in (None, UNKNOWN, [], ()):
            out(' '*6+'\\begin{EpydocFunctionRaises}\n')
            for name, descr in func_doc.exception_descrs:
                out(' '*10+'\\item[%s]\n\n' %
                    plaintext_to_latex('%s' % name))
                out(self.docstring_to_latex(descr, 10))
            out(' '*6+'\\end{EpydocFunctionRaises}\n\n')

        ## Overrides
        if var_doc.overrides not in (None, UNKNOWN):
            out('\\EpydocFunctionOverrides')
            if (func_doc.docstring in (None, UNKNOWN) and
                var_doc.overrides.value.docstring not in (None, UNKNOWN)):
                out('[1]')
            out('{%s}\n\n' 
                % self._hyperlink(var_doc.overrides, 
                                  '%s' % var_doc.overrides.canonical_name))

        # Add version, author, warnings, requirements, notes, etc.
        self.write_standard_fields(out, func_doc)

        out('    \\end{EpydocFunction}\n\n')

    def function_signature(self, var_doc):
        func_doc = var_doc.value
        func_name = var_doc.name
        
        s = '{%s}{' % self._hypertarget(var_doc, func_name)
        
        # This should never happen, but just in case:
        if func_doc not in (None, UNKNOWN):
            if func_doc.posargs == UNKNOWN:
                args = ['\\GenericArg{}']
            else:
                args = [self.func_arg(name, default) for (name, default)
                        in zip(func_doc.posargs, func_doc.posarg_defaults)]
            if func_doc.vararg:
                if func_doc.vararg == '...':
                    args.append('\\GenericArg{}')
                else:
                    args.append('\\VarArg{%s}' % plaintext_to_latex(func_doc.vararg))
            if func_doc.kwarg:
                args.append('\\KWArg{%s}' % plaintext_to_latex(func_doc.kwarg))
            
        s += '%\n    \\and '.join(args)
        s += '}%\n'
        
        return s

    def func_arg(self, name, default):
        s = '\\Param'
        
        if default is not None:
            s += "[%s]" % default.summary_pyval_repr().to_latex(None)
        s += '{%s}' % self._arg_name(name)

        return s

    def _arg_name(self, arg):
        if isinstance(arg, basestring):
            return plaintext_to_latex(arg)
        else:
            return '\\ArgList{%s}' % '\\and '.join([self._arg_name(a) for a in arg])

    #////////////////////////////////////////////////////////////
    #{ Variable List
    #////////////////////////////////////////////////////////////

    # Also used for the property list.
    def write_var_list(self, out, heading, doc, value_type, seclevel=1):
        groups = [(plaintext_to_latex(group_name),
                   doc.select_variables(group=group_name, imported=False,
                                        value_type=value_type,
                                        public=self._public_filter))
                  for group_name in doc.group_names()]

        # Discard any empty groups; and return if they're all empty.
        groups = [(g,vars) for (g,vars) in groups if vars]
        if not groups: return

        # Write a header.
        self.write_start_of(out, heading)
        out('  '+self.section(heading, seclevel))

        out('\\begin{Epydoc%sList}\n' % ''.join(heading.split()))

        for name, var_docs in groups:
            if name:
                out('\\EpydocInternalHeader{%s}\n' % name)
            for var_doc in var_docs:
                if isinstance(var_doc, PropertyDoc):
                    self.write_property_list_line(out, var_doc)
                else:
                    self.write_var_list_line(out, var_doc)
            # [xx] deal with inherited methods better????
            #if (self._inheritance == 'listed' and
            #    isinstance(container, ClassDoc)):
            #    out(self._inheritance_list(group, container.uid()))

        out('\\end{Epydoc%sList}\n\n' % ''.join(heading.split()))
    
    def write_var_list_line(self, out, var_doc):
        has_descr = var_doc.descr not in (None, UNKNOWN)
        has_type = var_doc.type_descr not in (None, UNKNOWN)
        has_repr = (var_doc.value not in (None, UNKNOWN) and
                    (var_doc.value.parse_repr is not UNKNOWN or
                     var_doc.value.pyval_repr() is not UNKNOWN))
                     
        out('\\EpydocVariable{%s}{' % self._hypertarget(var_doc, 
                                                        var_doc.name))
        if has_descr:
            out(self.docstring_to_latex(var_doc.descr, 10).strip())
        out('}{')
        if has_repr:
            out(var_doc.value.summary_pyval_repr().to_latex(None))
        out('}{')
        if has_type:
            out(self.docstring_to_latex(var_doc.type_descr, 12).strip())
        out('}\n')

    def write_property_list_line(self, out, var_doc):
        prop_doc = var_doc.value
        out('\\EpydocProperty{%s}{%s}{' 
            % (self._hypertarget(var_doc, var_doc.shortname), 
               plaintext_to_latex(var_doc.name, nbsp=True, breakany=True)))

        
        has_descr = prop_doc.descr not in (None, UNKNOWN)
        has_type = prop_doc.type_descr not in (None, UNKNOWN)
        if has_descr:
            out(self.docstring_to_latex(prop_doc.descr, 10).strip())
        out('}{')
        if has_type:
            ptype = self.docstring_to_latex(prop_doc.type_descr, 12).strip()
        # [xx] List the fget/fset/fdel functions?
        out('}\n')

    #////////////////////////////////////////////////////////////
    #{ Standard Fields
    #////////////////////////////////////////////////////////////

    # Copied from HTMLWriter:
    def write_standard_fields(self, out, doc):
        fields = []
        field_values = {}
        
        #if _sort_fields: fields = STANDARD_FIELD_NAMES [XX]
        
        for (field, arg, descr) in doc.metadata:
            if field not in field_values:
                fields.append(field)
            if field.takes_arg:
                subfields = field_values.setdefault(field,{})
                subfields.setdefault(arg,[]).append(descr)
            else:
                field_values.setdefault(field,[]).append(descr)

        for field in fields:
            if field.takes_arg:
                for arg, descrs in field_values[field].items():
                    self.write_standard_field(out, doc, field, descrs, arg)
                                              
            else:
                self.write_standard_field(out, doc, field, field_values[field])

    def write_standard_field(self, out, doc, field, descrs, arg=''):
        singular = field.singular
        plural = field.plural
        if arg:
            singular += ' (%s)' % arg
            plural += ' (%s)' % arg
        out(self._descrlist([self.docstring_to_latex(d) for d in descrs],
                            field.singular, field.plural, field.short))
            
    def _descrlist(self, items, singular, plural=None, short=0):
        if plural is None: plural = singular
        if len(items) == 0: return ''
        if len(items) == 1 and singular is not None:
            return '\\EpydocDescriptionItem{%s}{%s}\n\n' % (singular, items[0])
        if short:
            s = '\\begin{EpydocDescriptionShortList}{%s}%%\n    ' % plural
            s += '%\n    \\and '.join([item.strip() for item in items])
            s += '%\n\\end{EpydocDescriptionShortList}\n\n'
            return s
        else:
            s = '\\begin{EpydocDescriptionLongList}{%s}%%\n' % plural
            s += '\n\n'.join(['  \item %s%%' % item for item in items])
            s += '\n\\end{EpydocDescriptionLongList}\n\n'
            return s


    #////////////////////////////////////////////////////////////
    #{ Docstring -> LaTeX Conversion
    #////////////////////////////////////////////////////////////

    # We only need one linker, since we don't use context:
    class _LatexDocstringLinker(markup.DocstringLinker):
        def translate_indexterm(self, indexterm):
            indexstr = re.sub(r'["!|@]', r'"\1', indexterm.to_latex(self))
            return ('\\index{%s}\\textit{%s}' % (indexstr, indexstr))
        def translate_identifier_xref(self, identifier, label=None):
            if label is None: label = markup.plaintext_to_latex(identifier)
            return '\\texttt{%s}' % label
    _docstring_linker = _LatexDocstringLinker()
    
    def docstring_to_latex(self, docstring, indent=0, breakany=0):
        if docstring is None: return ''
        return docstring.to_latex(self._docstring_linker, indent=indent,
                                  hyperref=self._hyperref)
    
    #////////////////////////////////////////////////////////////
    #{ Helpers
    #////////////////////////////////////////////////////////////

    def write_header(self, out, where):
        out('%\n% API Documentation')
        if self._prj_name: out(' for %s' % self._prj_name)
        if isinstance(where, APIDoc):
            out('\n%% %s %s' % (self.doc_kind(where), where.canonical_name))
        else:
            out('\n%% %s' % where)
        out('\n%%\n%% Generated by epydoc %s\n' % epydoc.__version__)
        out('%% [%s]\n%%\n' % time.asctime(time.localtime(time.time())))

    def write_start_of(self, out, section_name):
        out('\n' + 75*'%' + '\n')
        out('%%' + section_name.center(71) + '%%\n')
        out(75*'%' + '\n\n')

    def section(self, title, depth=0, ref=None):
        sec = self.SECTIONS[depth+self._top_section]
        text = (('%s\n\n' % sec) % title)
        if ref:
            text += self._hypertarget(ref, "")
        return text
    
    def sectionstar(self, title, depth, ref=None):
        sec = self.STARSECTIONS[depth+self._top_section]
        text = (('%s\n\n' % sec) % plaintext_to_latex(title))
        if ref:
            text += self._hypertarget(ref, "")
        return text

    def doc_kind(self, doc):
        if isinstance(doc, ModuleDoc) and doc.is_package == True:
            return 'Package'
        elif (isinstance(doc, ModuleDoc) and
              doc.canonical_name[0].startswith('script')):
            return 'Script'
        elif isinstance(doc, ModuleDoc):
            return 'Module'
        elif isinstance(doc, ClassDoc):
            return 'Class'
        elif isinstance(doc, ClassMethodDoc):
            return 'Class Method'
        elif isinstance(doc, StaticMethodDoc):
            return 'Static Method'
        elif isinstance(doc, RoutineDoc):
            if isinstance(self.docindex.container(doc), ClassDoc):
                return 'Method'
            else:
                return 'Function'
        else:
            return 'Variable'

    def indexterm(self, doc, pos='only'):
        """Mark a term or section for inclusion in the index."""
        if not self._index: return ''
        if isinstance(doc, RoutineDoc) and not self._index_functions:
            return ''

        pieces = []
        
        if isinstance(doc, ClassDoc):
            classCrossRef = '\\index{\\EpydocIndex[%s]{%s}|see{%%s}}\n' \
              % (self.doc_kind(doc).lower(), 
                 self._dotted('%s' % doc.canonical_name))
        else:
            classCrossRef = None

        while isinstance(doc, ClassDoc) or isinstance(doc, RoutineDoc):
            if doc.canonical_name == UNKNOWN:
                return '' # Give up.
            pieces.append('\\EpydocIndex[%s]{%s}' %
                          (self.doc_kind(doc).lower(), 
                           self._dotted('%s' % doc.canonical_name)))
            doc = self.docindex.container(doc)
            if doc == UNKNOWN:
                return '' # Give up.

        pieces.reverse()
        if pos == 'only':
            term = '\\index{%s}\n' % '!'.join(pieces)
        elif pos == 'start':
            term = '\\index{%s|(}\n' % '!'.join(pieces)
        elif pos == 'end':
            term = '\\index{%s|)}\n' % '!'.join(pieces)
        else:
            raise AssertionError('Bad index position %s' % pos)
        
        if pos in ['only', 'start'] and classCrossRef is not None:
            term += classCrossRef % ('\\EpydocIndex[%s]{%s}' % (self.doc_kind(doc).lower(), 
                                                                self._dotted('%s'%doc.canonical_name)))
            
        return term

    def label(self, doc):
        return ':'.join(doc.canonical_name)
        
    def _hyperlink(self, target, name):
        return '\\EpydocHyperlink{%s}{%s}' % (self.label(target), self._dotted(name))

    def _hypertarget(self, uid, sig):
        return '\\EpydocHypertarget{%s}{%s}' % (self.label(uid), self._dotted(sig))

    def _dotted(self, name):
        return '\\EpydocDottedName{%s}' % name




    #: Map the Python encoding representation into mismatching LaTeX ones.
    latex_encodings = {
        'utf-8': 'utf8',
    }

    def get_latex_encoding(self):
        """
        @return: The LaTeX representation of the selected encoding.
        @rtype: C{str}
        """
        enc = self._encoding.lower()
        return self.latex_encodings.get(enc, enc)
