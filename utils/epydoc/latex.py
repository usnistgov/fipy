#!/usr/bin/python2.2
#
# epydoc.py: latex output
# Edward Loper
#
# Created [01/30/01 05:18 PM]
# $Id$
#

"""
Documentation to LaTeX converter.  This module defines a single class,
L{LatexFormatter}, which translates the API documentation encoded in a
L{DocMap} into a set of LaTeX files.

@var _LATEX_HEADER: The header for standard documentation LaTeX pages.
@var _HRULE: LaTeX code for a horizontal rule across the page
@var _SECTIONS: A list of string patterns that encode each numbered
    section level.
@var _STARSECTIONS: A list of string patterns that encode each
    un-numbered section level.

@bug: If a longtable is generated with a single row that's larger than
    a page, then LaTeX is unable to generate the correct output.  This
    should only be a problem if a variable or property has a very long
    description.
"""
__docformat__ = 'epytext en'

##################################################
## Imports
##################################################

# system imports
import sys, xml.dom.minidom, os.path, time, types, re
import pprint

# epydoc imports
import epydoc
import epydoc.markup as markup
from epydoc.uid import UID, Link, findUID, make_uid
from epydoc.imports import import_module
from epydoc.objdoc import DocMap, ModuleDoc, FuncDoc
from epydoc.objdoc import ClassDoc, Var, Raise, ObjDoc

##################################################
## CONSTANTS
##################################################

# Packages:
#   - fullpage: Bigger margins (c.f. normal latex article style)
#   - alltt: a verbatim-like environment
#   - parskip: put space between paragraphs
#   - fancyheadings: put section names in headings
#   - boxedminipage: boxes around functions/methods
#   - makeidx: generate an index
#   - multirow: multirow cells in tabulars
#   - longtable: multi-page tables (for var lists)
#   - tocbibind: add the index to the table of contents
#   - amssymb: extra math symbols.
_LATEX_HEADER = r"""
\documentclass{article}
\usepackage{alltt, parskip, fancyheadings, boxedminipage}
\usepackage{makeidx, multirow, longtable, tocbibind, amssymb}
\usepackage{fullpage}
%\usepackage[headings]{fullpage}
\begin{document}

\setlength{\parindent}{0ex}
\setlength{\fboxrule}{2\fboxrule}
\newlength{\BCL} % base class length, for base trees.

\pagestyle{fancy}
\renewcommand{\sectionmark}[1]{\markboth{#1}{}}
\renewcommand{\subsectionmark}[1]{\markright{#1}}

\newenvironment{Ventry}[1]%
  {\begin{list}{}{%
    \renewcommand{\makelabel}[1]{\texttt{##1:}\hfil}%
    \settowidth{\labelwidth}{\texttt{#1:}}%
    \setlength{\leftmargin}{\labelsep}%
    \addtolength{\leftmargin}{\labelwidth}}}%
  {\end{list}}
""".strip()

_HRULE = '\\rule{\\textwidth}{0.5\\fboxrule}\n\n'

_SECTIONS = ['\\part{%s}', '\\chapter{%s}', '\\section{%s}',
             '\\subsection{%s}', '\\subsubsection{%s}',
             '\\textbf{%s}']
_STARSECTIONS = ['\\part*{%s}', '\\chapter*{%s}', '\\section*{%s}',
                 '\\subsection*{%s}', '\\subsubsection*{%s}',
                 '\\textbf{%s}']

##################################################
## Docstring Linking (Crossreferences)
##################################################

class _LatexDocstringLinker(markup.DocstringLinker):
    def translate_indexterm(self, indexterm):
        indexstr = re.sub(r'["!|@]', r'"\1', indexterm.to_latex(self))
        return ('\\index{%s}\\textit{%s}' % (indexstr, indexstr))
    def translate_identifier_xref(self, identifier, label=None):
        if label is None: label = markup.plaintext_to_latex(identifier)
        return '\\texttt{%s}' % label

##################################################
## Documentation -> Latex Conversion
##################################################

class LatexFormatter:
    """
    Documentation to LaTeX converter.  The API documentation produced
    by C{LatexFormatter} consists of a single LaTeX document, divided
    into several different files.  In particular, C{LatexFormatter}
    generates the following files:
    
      - X{api.tex}: The top-level LaTeX file.  This file imports the
        other files, to create a single unified document.  This is the
        file that you should run C{latex} on.
      - X{I{module}-module.tex}: The API documentation for a module.
        I{module} is the complete dotted name of the module, such as
        C{sys} or C{epydoc.epytext}.
      - X{I{class}-class.tex}: The API documentation for a class,
        exception, or type.  I{class} is the complete dotted name of
        the class, such as C{epydoc.epytext.Token} or C{array.ArrayType}.
        These class documentation files are only created if the
        C{list_classes_separately} option is used; otherwise, the
        documentation for each class is included in its module's
        documentation file.

    The methods C{write_module} and C{write_class} used to generate
    individual module and class documentation LaTeX files.  These
    files can then be included as chapters or sections of other LaTeX
    documents (with C{"\\include"}).  When using these methods, you
    may wish to disable the C{crossref} option, which will turn off
    crossreferencing betweeen modules and classes, since some of these
    crossreference links will be broken if you only include some of
    the API documentation as chapters or sections of your document.

    @ivar _docmap: The documentation map, encoding the objects that
        should be documented.
    @ivar _show_private: Whether to include show private objects in
    the documentation.
    """
    
    def __init__(self, docmap, **kwargs):
        """
        Construct a new LaTeX formatter, using the given documentation map.
        @param docmap: The documentation to output.
        @type docmap: L{DocMap}
        @param kwargs: Keyword arguments:
            - C{prj_name}: The name of the project.  Defaults to
              none.  (type=C{string})
            - C{private}: Whether to create documentation for private
              objects.  By default, private objects are documented.
              (type=C{boolean})
            - C{crossref}: Whether to create crossreference links
              between classes and modules.  By default, crossreference
              links are created.  (type=C{boolean})
            - C{index}: Whether to generate an index.  If you generate
              an index, you will need to run C{makeindex} to make the
              C{.idx} file.  By default, an index is generated.
              (type=C{boolean})
            - C{list_classes_separately}: Whether to list classes in
              separate chapters, or to include them as sections of
              their modules' chapters.  By default, they are not listed
              separately.  (type=C{boolean})
            - C{exclude}: Whether to exclude inherited objects and
              imported objects that are not defined by any of the
              modules that are being documented.  By default, these
              objects are excluded.  (type=C{boolean})
            - C{alphabetical}: Whether to list modules in alphabetical
              order or in the order that they were specified in on
              the command line.  By default, modules are listed in
              alphabetical order.  (type=C{boolean})
        """
        self._docmap = docmap

        # Process keyword arguments
        self._show_private = kwargs.get('private', 0)
        self._prj_name = kwargs.get('prj_name', None)
        self._crossref = kwargs.get('crossref', 1)
        self._index = kwargs.get('index', 1)
        self._list_classes_separately=kwargs.get('list_classes_separately',0)
        self._inheritance = kwargs.get('inheritance', 'listed')
        self._exclude = kwargs.get('exclude', 1)
        self._top_section = 2
        self._index_functions = 1
        self._hyperref = 1
        self._module_uids = kwargs.get('modules', [])[:]
        if kwargs.get('alphabetical', 1):
            self._module_uids.sort()

    def write(self, directory=None, progress_callback=None):
        """
        Write the API documentation for the entire project to the
        given directory.

        @type directory: C{string}
        @param directory: The directory to which output should be
            written.  If no directory is specified, output will be
            written to the current directory.  If the directory does
            not exist, it will be created.
        @type progress_callback: C{function}
        @param progress_callback: A callback function that is called
            before each file is written, with the name of the created
            file.
        @rtype: C{None}
        @raise OSError: If C{directory} cannot be created,
        @raise OSError: If any file cannot be created or written to.
        """
        if not directory: directory = os.curdir
        
        # Create dest directory, if necessary
        if not os.path.isdir(directory):
            if os.path.exists(directory):
                raise OSError('%r is not a directory' % directory)
            os.mkdir(directory)

        # Write the module & class files.
        objects = self._filter(self._docmap.keys())
        objects.sort()
        for uid in objects:
            if uid.is_module():
                filename = os.path.join(directory, ('%s-module.tex' %
                                                    uid.name()))
                if progress_callback: progress_callback(filename)
                open(filename, 'w').write(self._module_to_latex(uid))
            elif uid.is_class() and self._list_classes_separately:
                if self._excluded(uid): continue
                filename = os.path.join(directory, ('%s-class.tex' %
                                                    uid.name()))
                if progress_callback: progress_callback(filename)
                open(filename, 'w').write(self._class_to_latex(uid))

        # Write the top-level file.
        filename = os.path.join(directory, 'api.tex')
        if progress_callback: progress_callback(filename)
        open(filename, 'w').write(self._topfile())

    def write_module(self, uid, filename):
        """
        Write the API documentation for the given module to
        C{filename}.
        @param uid: The unique identifier of the module to document.
        @type uid: L{UID}
        @param filename: The name of the file to write the
            documentation to.
        @type filename: C{string}
        @raise OSError: If C{directory} cannot be created,
        @raise OSError: If any file cannot be created or written to.
        @raise ValueError: If C{uid} is not the identifier for a module.
        """
        if not uid.is_module():
            raise ValueError('%s is not a module' % uid)
        open(filename, 'w').write(self._module_to_latex(uid))

    def write_class(self, uid, filename):
        """
        Write the API documentation for the given class to
        C{filename}.
        @param uid: The unique identifier of the class to document.
        @type uid: L{UID}
        @param filename: The name of the file to write the
            documentation to.
        @type filename: C{string}
        @raise OSError: If C{directory} cannot be created,
        @raise OSError: If any file cannot be created or written to.
        @raise ValueError: If C{uid} is not the identifier for a class.
        """
        if not uid.is_class():
            raise ValueError('%s is not a class' % uid)
        open(filename, 'w').write(self._class_to_latex(uid))

    def num_files(self):
        """
        @return: The number of files that this C{LatexFormatter} will
            generate.
        @rtype: C{int}
        """
        n = 1
        for uid in self._docmap.keys():
            if uid.is_private() and not self._show_private: continue
            if uid.is_module(): n += 1
            elif uid.is_class() and self._list_classes_separately:
                if not self._excluded(uid): n += 1
        return n
        
    #////////////////////////////////////////////////////////////
    # Main Doc File
    #////////////////////////////////////////////////////////////

    def _topfile(self):
        str = self._header('Inclue File')

        str += self._start_of('Header')
        str += _LATEX_HEADER + '\n'

        if self._index:
            str = re.sub(r'(\\begin{document})', '\\makeindex\n\\1', str)

        if self._hyperref:
            hyperref = (r'\\usepackage[usenames]{color}\n' +
                        r'\\definecolor{darkblue}{rgb}{0,0.05,0.35}\n' +
                        r'\\usepackage[dvips, pagebackref, ' +
                        'pdftitle={%s}, ' % (self._prj_name or '') +
                        'pdfcreator={epydoc %s}, ' % epydoc.__version__ +
                        'bookmarks=true, bookmarksopen=false, '+
                        'pdfpagemode=UseOutlines, colorlinks=true, '+
                        'linkcolor=black, anchorcolor=black, '+
                        'citecolor=black, filecolor=black, '+
                        'menucolor=black, pagecolor=black, '+
                        'urlcolor=darkblue]{hyperref}\n')
            str = re.sub(r'(\\begin{document})',
                         hyperref + '\\1', str)

        str += self._start_of('Title')
        str += '\\title{%s}\n' % self._text_to_latex(self._prj_name, 1)
        str += '\\author{API Documentation}\n'
        str += '\\maketitle\n'
        
        str += self._start_of('Table of Contents')
        str += '\\addtolength{\\parskip}{-1ex}\n'
        str += '\\tableofcontents\n'
        str += '\\addtolength{\\parskip}{1ex}\n'

        str += self._start_of('Includes')
        for uid in self._module_uids:
            str += '\\include{%s-module}\n' % uid.name()

        # If we're listing classes separately, put them after all the
        # modules.
        if self._list_classes_separately:
            uids = self._filter(self._docmap.keys())
            uids.sort()
            for uid in uids:
                if uid.is_class():
                    if self._excluded(uid): continue
                    str += '\\include{%s-class}\n' % uid.name()

        str += self._start_of('Index')
        str += '\\printindex\n\n'
        str += self._start_of('Footer')
        str += '\\end{document}\n\n'
        return str

    #////////////////////////////////////////////////////////////
    # Chapters
    #////////////////////////////////////////////////////////////

    def _module_to_latex(self, uid):
        # Get the module's documentation.
        doc = self._docmap[uid]

        # Start the chapter.
        str = self._header(uid)
        str += self._start_of('Module Description')
        str += '    ' + self._indexterm(uid, 'start')
        if uid.is_package():
            str += self._section('Package %s' % uid.name(), 0)
        else:
            str += self._section('Module %s' % uid.name(), 0)
        str += '    \\label{%s}\n' % self._uid_to_label(uid)

        # The module's description.
        if doc.descr():
            str += self._docstring_to_latex(doc.descr())

        # Add version, author, warnings, requirements, notes, etc.
        str += self._standard_fields(doc)

        # If it's a package, list the sub-modules.
        if doc.ispackage() and doc.modules():
            str += self._module_list(doc, doc.modules())

        # Class list. !! add summaries !!
        if self._list_classes_separately and doc.classes():
            str += self._class_list(doc, doc.classes())
                
        # Function List
        str += self._func_list(doc, doc.functions())

        # Variable list.
        if doc.variables():
            str += self._var_list(doc, doc.variables())

        # Class list.
        if not self._list_classes_separately:
            classes = self._filter(doc.classes())
            for cls in classes:
                str += self._class_to_latex(cls.target())
        
        str += '    ' + self._indexterm(uid, 'end')
        return str

    def _class_to_latex(self, uid):
        # Get the module's documentation.
        doc = self._docmap[uid]

        # Start the chapter.
        str = ''
        if self._list_classes_separately: str += self._header(uid)
        str += '    ' + self._indexterm(uid, 'start')
        str += self._start_of('Class Description')
        if self._list_classes_separately:
            seclevel = 0
            str += self._section('Class %s' % uid.name(), seclevel)
        else:
            seclevel = 1
            str += self._section('Class %s' % uid.shortname(), seclevel)
        str += '    \\label{%s}\n' % self._uid_to_label(uid)

        # The class base tree.
        if doc.bases():
            str += self._base_tree(uid)

        # The class's known subclasses.
        if doc.subclasses():
            str += self._subclasses(doc.subclasses(), uid)

        # The class's description
        if doc.descr():
            str += self._docstring_to_latex(doc.descr())
        
        # Add version, author, warnings, requirements, notes, etc.
        str += self._standard_fields(doc)

        # Methods.
        str += self._func_list(doc, doc.methods(),
                               'Methods', seclevel+1)
        str += self._func_list(doc, doc.staticmethods(),
                               'Static Methods', seclevel+1)
        str += self._func_list(doc, doc.classmethods(),
                               'Class Methods', seclevel+1)

        if doc.properties():
            str += self._property_list(doc, doc.properties(), 
                                       'Properties', seclevel+1)
        if doc.ivariables():
            str += self._var_list(doc, doc.ivariables(), 
                                  'Instance Variables', seclevel+1)
        if doc.cvariables():
            str += self._var_list(doc, doc.cvariables(), 
                                  'Class Variables', seclevel+1)

        # End mark for the class's index entry.
        str += '    ' + self._indexterm(uid, 'end')
        
        return str

    #////////////////////////////////////////////////////////////
    # Class List
    #////////////////////////////////////////////////////////////
    def _class_list(self, container, classes):
        classes = self._filter(classes)
        if len(classes) == 0: return ''
        
        groups = container.by_group(classes)

        str = self._start_of('Classes')
        str += self._section('Classes', 1)
        str += '\\begin{itemize}'
        str += '  \\setlength{\\parskip}{0ex}\n'
        
        # Create the portion of the table containing the group
        # entries.  Do this first, so we can see what's not in any
        # group; but add it to the string last, so the groupless
        # properties are at the top.
        for name, group in groups:
            if name is not None:
                str += '  \\item \\textbf{%s}\n' % name
                str += '  \\begin{itemize}\n'
            # Add the lines for each class
            for cls in group:
                str += self._class_list_line(cls)
            if name is not None:
                str += '  \end{itemize}\n'
        
        return str + '\\end{itemize}'

    def _class_list_line(self, link):
        cname = link.name()
        cls = link.target()
        if not self._docmap.has_key(cls): return ''
        cdoc = self._docmap[cls]
        str = '  ' + '\\item \\textbf{'
        str += self._text_to_latex(cname) + '}'
        if cdoc and cdoc.descr():
            str += ': %s\n' % self._summary(cdoc, cls.module())
        if self._crossref:
            str += ('\n  \\textit{(Section \\ref{%s}' %
                    self._uid_to_label(cls))
            str += (', p.~\\pageref{%s})}\n\n' %
                    self._uid_to_label(cls))
        return str
        
    #////////////////////////////////////////////////////////////
    # Property List
    #////////////////////////////////////////////////////////////
    def _property_list(self, container, properties, 
                       heading='Properties', seclevel=1):
        properties = self._filter(properties)
        if len(properties) == 0: return ''

        groups = container.by_group(properties)
        str = self._start_of(heading)
        str += '  '+self._section(heading, seclevel)

        str += '\\begin{longtable}'
        str += '{|p{.30\\textwidth}|'
        str += 'p{.62\\textwidth}|l}\n'
        str += '\\cline{1-2}\n'

        # Set up the headers & footer (this makes the table span
        # multiple pages in a happy way).
        str += '\\cline{1-2} '
        str += '\\centering \\textbf{Name} & '
        str += '\\centering \\textbf{Description}& \\\\\n'
        str += '\\cline{1-2}\n'
        str += '\\endhead'
        str += '\\cline{1-2}'
        str += '\\multicolumn{3}{r}{\\small\\textit{'
        str += 'continued on next page}}\\\\'
        str += '\\endfoot'
        str += '\\cline{1-2}\n'
        str += '\\endlastfoot'

        # Create the portion of the table containing the group
        # entries.  Do this first, so we can see what's not in any
        # group; but add it to the string last, so the groupless
        # properties are at the top.
        for name, group in groups:
            # Print a header within the table
            if name is not None:
                str += '\\multicolumn{2}{|l|}{'
                str += '\\textbf{%s}}\\\\\n' % name
                str += '\\cline{1-2}\n'
            # Add the lines for each property
            for property in group:
                str += self._property_list_line(property, container.uid())
            if self._inheritance == 'listed':
                str += self._inheritance_list_line(group, container.uid())
        return str + '\\end{longtable}\n\n'

    def _property_list_line(self, link, container):
        inherit = (container.is_class() and container != link.target().cls())
        if inherit and self._inheritance == 'listed': return ''
            
        prop = link.target()
        pdoc = self._docmap.get(prop)
        if pdoc is None: return ''
        str = '\\raggedright '
        str += self._text_to_latex(prop.shortname(), 1, 1)
        str += ' & '

        if pdoc.descr():
            str += '\\raggedright '
            str += self._docstring_to_latex(pdoc.descr(), 10).strip()
        str += '&\\\\\n'
        str += '\\cline{1-2}\n'
        return str

    #////////////////////////////////////////////////////////////
    # Variable List
    #////////////////////////////////////////////////////////////

    def _var_list(self, container, variables, 
                  heading='Variables', seclevel=1):
        variables = self._filter(variables)
        if len(variables) == 0: return ''
        
        groups = container.by_group(variables)

        str = self._start_of(heading)
        str += '  '+self._section(heading, seclevel)

        str += '\\begin{longtable}'
        str += '{|p{.30\\textwidth}|'
        str += 'p{.62\\textwidth}|l}\n'
        str += '\\cline{1-2}\n'

        # Set up the headers & footer (this makes the table span
        # multiple pages in a happy way).
        str += '\\cline{1-2} '
        str += '\\centering \\textbf{Name} & '
        str += '\\centering \\textbf{Description}& \\\\\n'
        str += '\\cline{1-2}\n'
        str += '\\endhead'
        str += '\\cline{1-2}'
        str += '\\multicolumn{3}{r}{\\small\\textit{'
        str += 'continued on next page}}\\\\'
        str += '\\endfoot'
        str += '\\cline{1-2}\n'
        str += '\\endlastfoot'

        # Create the portion of the table containing the group
        # entries.  Do this first, so we can see what's not in any
        # group; but add it to the string last, so the groupless
        # functions are at the top.
        for name, group in groups:
            # Print a header within the table
            if name is not None:
                str += '\\multicolumn{2}{|l|}{'
                str += '\\textbf{%s}}\\\\\n' % name
                str += '\\cline{1-2}\n'
            # Add the lines for each variable
            for var in group:
                str += self._var_list_line(var, container.uid())
            if (self._inheritance == 'listed' and
                isinstance(container, ClassDoc)):
                str += self._inheritance_list_line(group, container.uid())
        return str + '\\end{longtable}\n\n'
    
    def _var_list_line(self, var, container):
        inherit = (container.is_class() and container != var.uid().cls())
        if inherit and self._inheritance == 'listed': return ''
            
        str = '\\raggedright '
        str += self._text_to_latex(var.name(), 1, 1) + ' & '

        if var.descr() or var.has_value():
            str += '\\raggedright '
        if var.descr():
            str += self._docstring_to_latex(var.descr(), 10).strip()
            if var.has_value() or var.type(): str += '\n\n'
        if var.has_value():
            str += '\\textbf{Value:} \n'
            str += self._pprint_var_value(var, 80)
        if var.type():
            ptype = self._docstring_to_latex(var.type(), 12).strip()
            str += '%s\\textit{(type=%s)}' % (' '*12, ptype)
        str += '&\\\\\n'
        str += '\\cline{1-2}\n'
        return str

    def _pprint_var_value(self, var, maxwidth=100):
        val = var.uid().value()
        try: val = `val`
        except: val = '...'
        if len(val) > maxwidth: val = val[:maxwidth-3] + '...'
        if '\n' in val:
            return ('\\begin{alltt}\n%s\\end{alltt}' %
                    self._text_to_latex(val, 0, 1))
        else:
            return '{\\tt %s}' % self._text_to_latex(val, 1, 1)
    
    #////////////////////////////////////////////////////////////
    # Function List
    #////////////////////////////////////////////////////////////
    
    def _func_list(self, container, functions, 
                   heading='Functions', seclevel=1):
        functions = self._filter(functions)
        if len(functions) == 0: return ''

        groups = container.by_group(functions)

        str = self._start_of(heading)
        str += '  '+self._section(heading, seclevel)

        # Create the portion of the table containing the group
        # entries.  Do this first, so we can see what's not in any
        # group; but add it to the string last, so the groupless
        # functions are at the top.
        for name, group in groups:
            # Print a header within the table
            if name is not None:
                str += '\n%s\\large{%s}\n' % (_HRULE, name)
            # Add the lines for each function
            for link in group:
                str += self._func_list_box(link, container.uid())
            if (self._inheritance == 'listed' and
                isinstance(container, ClassDoc)):
                str += self._inheritance_list(group, container.uid())

        return str
            
    def _func_list_box(self, link, cls):
        str = ''
        fname = link.name()
        fuid = link.target()
        if fuid.is_method() or fuid.is_builtin_method():
            container = fuid.cls()
            # (If container==ClassType, it's (probably) a class method.)
            inherit = (container != cls and
                       container.value() is not types.ClassType)
        else:
            inherit = 0
            try: container = fuid.module()
            except TypeError: container = None

        # Don't include inherited functions, if inheritance=listed.
        if inherit and self._inheritance == 'listed': return ''

        # If we don't have documentation for the function, then we
        # can't say anything about it.
        if not self._docmap.has_key(fuid): return ''
        fdoc = self._docmap[fuid]

        # What does this method override?
        foverrides = fdoc.overrides()

        # Try to find a documented ancestor.
        if fuid.is_any_method() and not fdoc.has_docstring():
            inhdoc = self._docmap.documented_ancestor(fuid) or fdoc
        else:
            inhdoc = fdoc
        inherit_docs = (inhdoc is not fdoc)

        # nb: this gives the containing section, not a reference
        # directly to the function.
        if not inherit:
            str += '    \\label{%s}\n' % self._uid_to_label(fuid)
        
        fsig = self._func_signature(fname, fdoc)
        if not inherit:
            str += '    ' + self._indexterm(fuid)
        str += '    \\vspace{0.5ex}\n\n'
        str += '    \\begin{boxedminipage}{\\textwidth}\n\n'
        str += '    %s\n\n' % fsig

        # Use the inherited docs for everything but the signature.
        fdoc = inhdoc

        if fdoc.has_docstring():
            str += '    \\vspace{-1.5ex}\n\n'
            str += '    \\rule{\\textwidth}{0.5\\fboxrule}\n'
        
        fdescr=fdoc.descr()
        fparam = fdoc.parameter_list()[:]
        freturn = fdoc.returns()
        fraises = fdoc.raises()
        
        # Don't list parameters that don't have any extra info.
        f = lambda p:p.descr() or p.type()
        fparam = filter(f, fparam)

        # Description
        if fdescr:
            str += self._docstring_to_latex(fdescr, 4)
            str += '    \\vspace{1ex}\n\n'

        # Parameters
        if fparam:
            longest = max([len(p.name()) for p in fparam])
            str += ' '*6+'\\textbf{Parameters}\n'
            str += ' '*6+'\\begin{quote}\n'
            str += '        \\begin{Ventry}{%s}\n\n' % (longest*'x')
            for param in fparam:
                if param.listed_under(): continue
                str += ' '*10+'\\item[' + self._text_to_latex(param.name())
                if param.shared_descr_params():
                    for p in param.shared_descr_params():
                        str += ', %s' % self._text_to_latex(p.name())
                str += ']\n\n'
                if param.descr():
                    str += self._docstring_to_latex(param.descr(), 10)
                if param.shared_descr_params():
                    for p in [param]+param.shared_descr_params():
                        if not p.type(): continue
                        ptype = self._docstring_to_latex(p.type(), 14).strip()
                        str += (' '*12+'\\textit{(typeof %s=%s)}\n\n' %
                                (self._text_to_latex(p.name()), ptype))
                elif param.type():
                    ptype = self._docstring_to_latex(param.type(), 12).strip()
                    str += ' '*12+'\\textit{(type=%s)}\n\n' % ptype
            str += '        \\end{Ventry}\n\n'
            str += ' '*6+'\\end{quote}\n\n'
            str += '    \\vspace{1ex}\n\n'

        # Returns
        if freturn.descr() or freturn.type():
            str += ' '*6+'\\textbf{Return Value}\n'
            str += ' '*6+'\\begin{quote}\n'
            if freturn.descr():
                str += self._docstring_to_latex(freturn.descr(), 6)
                if freturn.type():
                    rtype = self._docstring_to_latex(freturn.type(), 6).strip()
                    str += ' '*6+'\\textit{(type=%s)}\n\n' % rtype
            elif freturn.type():
                str += self._docstring_to_latex(freturn.type(), 6)
            str += ' '*6+'\\end{quote}\n\n'
            str += '    \\vspace{1ex}\n\n'

        # Raises
        if fraises:
            str += ' '*6+'\\textbf{Raises}\n'
            str += ' '*6+'\\begin{quote}\n'
            str += '        \\begin{description}\n\n'
            for fraise in fraises:
                str += '          '
                str += '\\item[\\texttt{'+fraise.name()+'}]\n\n'
                str += self._docstring_to_latex(fraise.descr(), 10)
            str += '        \\end{description}\n\n'
            str += ' '*6+'\\end{quote}\n\n'
            str += '    \\vspace{1ex}\n\n'

        ## Overrides
        if foverrides:
            str += ('      Overrides: %s' %
                    self._text_to_latex(foverrides.name()))
            if inherit_docs:
                str += ' \textit{(inherited documentation)}'
            str += '\n\n'

        # Add version, author, warnings, requirements, notes, etc.
        str += self._standard_fields(fdoc)

        str += '    \\end{boxedminipage}\n\n'
        return str

    def _func_signature(self, fname, fdoc, show_defaults=1):
        str = '\\raggedright '
        str += '\\textbf{%s}' % self._text_to_latex(fname)
        str += '('
        str += self._params_to_latex(fdoc.parameters(), show_defaults)
        
        if fdoc.vararg():
            vararg_name = self._text_to_latex(fdoc.vararg().name())
            vararg_name = '\\textit{%s}' % vararg_name
            if vararg_name != '\\textit{...}':
                vararg_name = '*%s' % vararg_name
            str += '%s, ' % vararg_name
        if fdoc.kwarg():
            str += ('**\\textit{%s}, ' %
                    self._text_to_latex(fdoc.kwarg().name()))
        if str[-1] != '(': str = str[:-2]

        return str + ')'
    
    def _params_to_latex(self, parameters, show_defaults):
        str = ''
        for param in parameters:
            if type(param) in (type([]), type(())):
                sublist = self._params_to_latex(param, show_defaults)
                str += '(%s), ' % sublist[:-2]
            else:
                str += '\\textit{%s}' % self._text_to_latex(param.name())
                if show_defaults and param.default() is not None:
                    default = param.default()
                    if len(default) > 60:
                        default = default[:57]+'...'
                    str += ('=\\texttt{%s}' %
                            self._text_to_latex(default, 1, 1))
                str += ', '
        return str

    #////////////////////////////////////////////////////////////
    # Inheritance lists
    #////////////////////////////////////////////////////////////
    
    def _inheritance_list(self, links, cls):
        # Group the objects by defining class
        inh_dict = {}
        for link in links:
            if isinstance(link, Link): key = link.target().cls()
            else: key = link.uid().cls()
            if key == cls: continue
            if key is None: continue
            if not inh_dict.has_key(key): inh_dict[key] = []
            inh_dict[key].append(link)

        if not inh_dict: return ''

        str = ''#'\\begin{boxedminipage}{\\textwidth}\n'
        inh_items = inh_dict.items()
        inh_items.sort(lambda a,b: cmp(a[0], b[0]))
        for (base, obj_links) in inh_items:
            str += ('  \\textbf{Inherited from %s:}\n' %
                    self._text_to_latex(base.shortname()))
            for link in obj_links:
                str += '    '
                str += self._text_to_latex(link.name())
                str += ',\n'
            str = str[:-2] + '\n    \\\\\n'
        return str[:-7] #+ '\\end{boxedminipage}\n'

    def _inheritance_list_line(self, links, cls):
        # Group the objects by defining class
        inh_dict = {}
        for link in links:
            if isinstance(link, Link): key = link.target().cls()
            else: key = link.uid().cls()
            if key == cls: continue
            if key is None: continue
            if not inh_dict.has_key(key): inh_dict[key] = []
            inh_dict[key].append(link)

        if not inh_dict: return ''

        str = ''
        inh_items = inh_dict.items()
        inh_items.sort(lambda a,b: cmp(a[0], b[0]))
        for (base, obj_links) in inh_items:
            str += '\\multicolumn{2}{|p{\\textwidth}|}{\n'
            str += ('  \\textbf{Inherited from %s:}\n' %
                    self._text_to_latex(base.shortname()))
            for link in obj_links:
                str += '    '
                str += self._text_to_latex(link.name())
                if self._crossref:
                    str += (' \\textit{(p.~\\pageref{%s})}' %
                            self._uid_to_label(base))
                str += ',\n'
            str = str[:-2] + '}\n    \\\\\n'
        return str + '\\cline{1-2}\n'

    #////////////////////////////////////////////////////////////
    # Docstring -> LaTeX Conversion
    #////////////////////////////////////////////////////////////

    _docstring_linker = _LatexDocstringLinker()
    def _docstring_to_latex(self, docstring, indent=0, breakany=0):
        if docstring is None: return ''
        return docstring.to_latex(self._docstring_linker, indent=indent,
                                  hyperref=self._hyperref)
    
    #////////////////////////////////////////////////////////////
    # Base class trees
    #////////////////////////////////////////////////////////////

    def _find_tree_width(self, uid):
        width = 2
        if self._docmap.has_key(uid):
            for base in self._docmap[uid].bases():
                width = max(width, self._find_tree_width(base.target())+2)

        return width

    def _base_tree(self, uid, width=None, linespec=None):
        if width is None:
            width = self._find_tree_width(uid)+2
            linespec = []
            str = ('&'*(width-4)+'\\multicolumn{2}{l}{\\textbf{%s}}\n' %
                   self._text_to_latex(uid.shortname()))
            str += '\\end{tabular}\n\n'
            top = 1
        else:
            str = self._base_tree_line(uid, width, linespec)
            top = 0
        
        bases = self._docmap[uid].bases()
        
        for i in range(len(bases)-1, -1, -1):
            base = bases[i].target()
            spec = (i > 0)
            str = self._base_tree(base, width, [spec]+linespec) + str

        if top:
            str = '\\begin{tabular}{%s}\n' % (width*'c') + str

        return str

    def _base_tree_line(self, uid, width, linespec):
        # linespec is a list of booleans.

        str = '%% Line for %s, linespec=%s\n' % (uid.name(), linespec)

        labelwidth = width-2*len(linespec)-2

        # The base class name.
        shortname = self._text_to_latex(uid.name())
        str += ('\\multicolumn{%s}{r}{' % labelwidth)
        str += '\\settowidth{\\BCL}{%s}' % shortname
        str += '\\multirow{2}{\\BCL}{%s}}\n' % shortname

        # The vertical bars for other base classes (top half)
        for vbar in linespec:
            if vbar: str += '&&\\multicolumn{1}{|c}{}\n'
            else: str += '&&\n'

        # The horizontal line.
        str += '  \\\\\\cline{%s-%s}\n' % (labelwidth+1, labelwidth+1)

        # The vertical bar for this base class.
        str += '  ' + '&'*labelwidth
        str += '\\multicolumn{1}{c|}{}\n'

        # The vertical bars for other base classes (bottom half)
        for vbar in linespec:
            if vbar: str += '&\\multicolumn{1}{|c}{}&\n'
            else: str += '&&\n'
        str += '  \\\\\n'

        return str
        
    #////////////////////////////////////////////////////////////
    # Module hierarchy trees
    #////////////////////////////////////////////////////////////
    
    def _module_tree_item(self, uid=None, depth=0):
        """
        Helper function for L{_module_tree} and L{_module_list}.
        
        @rtype: C{string}
        """
        if uid is None: return ''

        doc = self._docmap.get(uid, None)
        str = ' '*depth + '\\item \\textbf{'
        str += self._text_to_latex(uid.shortname()) +'}'
        if doc and doc.descr():
            str += ': %s\n' % self._summary(doc, uid)
        if self._crossref:
            str += ('\n  \\textit{(Section \\ref{%s}' %
                    self._uid_to_label(uid))
            str += ', p.~\\pageref{%s})}\n\n' % self._uid_to_label(uid)
        if doc and doc.ispackage() and doc.modules():
            str += ' '*depth + '  \\begin{itemize}\n'
            str += ' '*depth + '\\setlength{\\parskip}{0ex}\n'
            modules = [l.target() for l in self._filter(doc.modules())]
            for module in modules:
                str += self._module_tree_item(module, depth+4)
            str += ' '*depth + '  \\end{itemize}\n'
        return str

    def _module_tree(self, sortorder=None):
        """
        @return: The HTML code for the module hierarchy tree.  This is
            used by C{_trees_to_latex} to construct the hiearchy page.
            (Well, actually, it's not used by anything at present.)
        @rtype: C{string}
        """
        str = '\\begin{itemize}\n'
        str += '\\setlength{\\parskip}{0ex}\n'
        uids = self._filter(self._docmap.keys())
        uids.sort()
        #docs.sort(lambda a,b: cmp(a[0], b[0]))
        # Find all top-level packages. (what about top-level
        # modules?)
        for uid in uids:
            doc = self._docmap[uid]
            if not isinstance(doc, ModuleDoc): continue
            if not doc.package():
                str += self._module_tree_item(uid)
        return str +'\\end{itemize}\n'

    def _module_list(self, container, modules):
        """
        @return: The HTML code for the module hierarchy tree,
            containing the given modules.  This is used by
            L{_module_to_latex} to list the submodules of a package.
        @rtype: C{string}
        """
        modules = self._filter(modules)
        if len(modules) == 0: return ''
        str = self._start_of('Modules')
        str += self._section('Modules', 1)
        str += '\\begin{itemize}\n'
        str += '\\setlength{\\parskip}{0ex}\n'

        groups = container.by_group(modules)
        
        # Create the portion of the table containing the group
        # entries.  Do this first, so we can see what's not in any
        # group; but add it to the string last, so the groupless
        # properties are at the top.
        for name, group in groups:
            # Print a header within the table
            if name is not None:
                str += '  \\item \\textbf{%s}\n' % name
                str += '  \\begin{itemize}\n'
            # Add the lines for each module
            for link in group:
                str += self._module_tree_item(link.target())
            if name is not None:
                str += '  \end{itemize}\n'
        
        return str + '\\end{itemize}\n\n'

    #////////////////////////////////////////////////////////////
    # Helpers
    #////////////////////////////////////////////////////////////

    def _indexterm(self, uid, pos='only'):
        if not self._index: return ''
        if uid.is_routine() and not self._index_functions: return ''

        str = ''
        u = uid
        while (u.is_routine() or u.is_class()):
            str = '!%s \\textit{(%s)}%s' % (self._text_to_latex(u.shortname()),
                               self._kind(u).lower(), str)
            u = u.parent()

        str = '%s \\textit{(%s)}%s' % (self._text_to_latex(u.name()),
                          self._kind(u).lower(), str)

        if pos == 'only': return '\\index{%s}\n' % str
        elif pos == 'start': return '\\index{%s|(}\n' % str
        elif pos == 'end': return '\\index{%s|)}\n' % str
        else:
            raise AssertionError('Bad index position %s' % pos)

    def _text_to_latex(self, str, nbsp=0, breakany=0):
        """
        @param breakany: Insert hyphenation marks, so that LaTeX can
        break the resulting string at any point.  This is useful for
        small boxes (e.g., the type box in the variable list table).
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

    def _header(self, where):
        str = '%\n% API Documentation'
        if self._prj_name: str += ' for %s' % self._prj_name
        if isinstance(where, UID):
            str += '\n%% %s %s' % (self._kind(where), where.name())
        else:
            str += '\n%% %s' % where
        str += '\n%%\n%% Generated by epydoc %s\n' % epydoc.__version__
        str += '%% [%s]\n%%\n' % time.asctime(time.localtime(time.time()))
        return str

    def _kind(self, uid):
        if uid.is_package(): return 'Package'
        elif uid.is_module(): return 'Module'
        elif uid.is_class(): return 'Class'
        elif uid.is_method() or uid.is_builtin_method(): return 'Method'
        elif uid.is_routine(): return 'Function'
        elif uid.is_variable(): return 'Variable'
        else: raise AssertionError, 'Bad UID type for _name'

    def _section(self, title, depth):
        sec = _SECTIONS[depth+self._top_section]
        return (('%s\n\n' % sec) % self._text_to_latex(title))                
    
    def _sectionstar(self, title, depth):
        sec = _STARSECTIONS[depth+self._top_section]
        return (('%s\n\n' % sec) % self._text_to_latex(title))

    def _start_of(self, section_name):
        str = '\n' + 75*'%' + '\n'
        str += '%%' + ((71-len(section_name))/2)*' '
        str += section_name
        str += ((72-len(section_name))/2)*' ' + '%%\n'
        str += 75*'%' + '\n\n'
        return str

    def _uid_to_label(self, uid):
        return uid.name().replace('.', ':')
                
    def _cmp_name(self, name1, name2):
        """
        Compare uid1 and uid2 by their names, using the following rules: 
          - C{'__init__'} < anything.
          - public < private.
          - otherwise, sort alphabetically by name (ignoring case)
    
        @return: -1 if C{uid1<uid2}; 0 if C{uid1==uid2}; and 1 if
            C{uid1>uid2}.
        @rtype: C{int}
        """
        if (name2 == '__init__'): return 1
        if (name1 == '__init__'): return -1
        if name1 == name2: return 0
        if self._is_private(name1) and not self._is_private(name2): return 1
        if self._is_private(name2) and not self._is_private(name1): return -1
        return cmp(name1.lower(), name2.lower())
    
    def _is_private(self, str):
        """
        @return: true if C{str} is the name of a private Python object.
        @rtype: C{boolean}
        """
        if str == '...': return 0
        for piece in str.split('.'):
            if piece[:1] == '_' and piece[-1:] != '_': return 1
        return 0

    def _filter(self, links):
        """
        Filter a list of C{Link}s.  If L{_show_private} is false, then
        filter out all private objects; otherwise, perform no
        filtering.

        @param links: The list of C{Link}s to be filtered.
        @type links: C{list} of L{Link}
        @return: The filtered list of links.
        @rtype: C{list} of L{Link}
        """
        # Filter out private objects.
        if not self._show_private:
            return [l for l in links if l.is_public()]
        else:
            return links

    def _standard_fields(self, doc):
        """
        @return: HTML code containing descriptions of the epytext
        fields that are common to all L{ObjDoc}s (except for C{descr}).
        @rtype: C{string}
        @param doc: The object whose fields should be described.
        """
        uid = doc.uid()
        if uid.is_module() or uid.is_class(): container = uid
        else: container = uid.cls() or uid.module()
        str = ''

        for field in doc.fields():
            values = doc.field_values(field)
            if not values: continue
            items = [self._docstring_to_latex(v) for v in values]
            str += self._descrlist(items, field.singular,
                                   field.plural, field.short)
            
        return str
            
    def _descrlist(self, items, singular, plural=None, short=0):
        if plural is None: plural = singular
        if len(items) == 0: return ''
        if len(items) == 1 and singular is not None:
            return '\\textbf{%s:} %s\n\n' % (singular, items[0])
        if short:
            str = '\\textbf{%s:}\n' % plural
            items = [item.strip() for item in items]
            return str + ',\n    '.join(items) + '\n\n'
        else:
            str = '\\textbf{%s:}\n' % plural
            str += '\\begin{quote}\n'
            str += '  \\begin{itemize}\n\n  \item '
            str += '    \\setlength{\\parskip}{0.6ex}\n'
            str += '\n\n  \item '.join(items)
            return str + '\n\n\\end{itemize}\n\n\\end{quote}\n\n'

    def _subclasses(self, subclasses, container):
        """
        @return: The LaTeX code for the subclasses field.
        """
        items = [self._text_to_latex(sc.name()) for sc in subclasses]
        return self._descrlist(items, 'Known Subclasses', short=1)

    def _summary(self, doc, container=None):
        """
        @return: The LATEX code for the summary description of the
            object documented by C{doc}.  A summary description is the
            first sentence of the C{doc}'s 'description' field.  If the
            C{doc} has no 'description' field, but does have a
            'return' field, then the summary is taken from the return
            field instead.
        @rtype: C{string}
        @param doc: The documentation for the object whose summary
            should be returned.
        @type doc: L{objdoc.ObjDoc}
        @param container: The container object for C{doc}, or C{None}
            if there is none.  This container object is used to
            resolve links (E{L}{...}) in the epytext.
        @type container: L{uid.UID}
        """
        descr = doc.descr()

        # Try to find a documented ancestor.
        if isinstance(doc, FuncDoc):
            while (not doc.has_docstring() and doc.matches_override() and
                   self._docmap.has_key(doc.overrides())):
                doc = self._docmap[doc.overrides()]

        if descr != None:
            str = self._docstring_to_latex(descr.summary()).strip()
            return str
        elif (isinstance(doc, FuncDoc) and
              doc.returns().descr() is not None):
            summary = doc.returns().descr().summary()
            summary = self._docstring_to_latex(summary).strip()
            summary = summary[:1].lower() + summary[1:]
            return ('Return '+ summary)
        else:
            return ''

    def _excluded(self, x):
        """
        @return: True if the given object should be excluded from the
        documentation (since it was imported or inherited from a
        module that we're not documenting).
        """
        if not self._exclude: return 0
        if isinstance(x, Link): x = x.target()
        if isinstance(x, Var): x = x.var()
        if x.is_module(): return 0
        if x.module() is None: return 0
        return not self._docmap.has_key(x.module())
        
