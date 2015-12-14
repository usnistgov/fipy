import os.path
import posixpath

from docutils import io, nodes, statemachine, utils
from docutils.parsers.rst import directives, states

from sphinx import addnodes
from sphinx.util import patfilter, ws_re, caption_ref_re, url_re, docname_join
from sphinx.util.compat import Directive, directive_dwim, make_admonition

class FiPyExample(Directive):

    """
    Include content read from a separate source file.

    Content may be parsed by the parser, or included as a literal
    block.  The encoding of the included file can be specified.  Only
    a part of the given file argument may be included by specifying
    text to match before and/or after the text to be used.
    """

    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec = {'literal': directives.flag,
                   'encoding': directives.encoding,
                   'start-after': directives.unchanged_required,
                   'end-before': directives.unchanged_required}

    standard_include_path = os.path.join(os.path.dirname(states.__file__),
                                         'include')

    def run(self):
        """Include a reST file as part of the content of this reST file."""
        if not self.state.document.settings.file_insertion_enabled:
            raise self.warning('"%s" directive disabled.' % self.name)
        source = self.state_machine.input_lines.source(
            self.lineno - self.state_machine.input_offset - 1)
        source_dir = os.path.dirname(os.path.abspath(source))
        path = directives.path(self.arguments[0])
        if path.startswith('<') and path.endswith('>'):
            path = os.path.join(self.standard_include_path, path[1:-1])
        path = os.path.normpath(os.path.join(source_dir, path))
        path = utils.relative_path(None, path)
        encoding = self.options.get(
            'encoding', self.state.document.settings.input_encoding)
        try:
            self.state.document.settings.record_dependencies.add(path)
            include_file = io.FileInput(
                source_path=path, encoding=encoding,
                error_handler=(self.state.document.settings.\
                               input_encoding_error_handler),
                handle_io_errors=None)
        except IOError, error:
            raise self.severe('Problems with "%s" directive path:\n%s: %s.'
                              % (self.name, error.__class__.__name__, error))
        try:
            include_text = include_file.read()
        except UnicodeError, error:
            raise self.severe(
                'Problem with "%s" directive:\n%s: %s'
                % (self.name, error.__class__.__name__, error))
        # start-after/end-before: no restrictions on newlines in match-text,
        # and no restrictions on matching inside lines vs. line boundaries
        after_text = self.options.get('start-after', '"""')
        if after_text:
            # skip content in include_text before *and incl.* a matching text
            after_index = include_text.find(after_text)
            if after_index < 0:
                raise self.severe('Problem with "start-after" option of "%s" '
                                  'directive:\nText not found.' % self.name)
            include_text = include_text[after_index + len(after_text):]
        before_text = self.options.get('end-before', '"""')
        if before_text:
            # skip content in include_text after *and incl.* a matching text
            before_index = include_text.find(before_text)
            if before_index < 0:
                raise self.severe('Problem with "end-before" option of "%s" '
                                  'directive:\nText not found.' % self.name)
            include_text = include_text[:before_index]
        if 'literal' in self.options:
            literal_block = nodes.literal_block(include_text, include_text,
                                                source=path)
            literal_block.line = 1
            return [literal_block]
        else:
            include_lines = statemachine.string2lines(include_text,
                                                      convert_whitespace=1)

            self.state_machine.insert_input(include_lines, path)
            self.state_machine.insert_input("-" * len(path), path)
            self.state_machine.insert_input(path, path)
            return []

def exampletree_directive(name, arguments, options, content, lineno,
                          content_offset, block_text, state, state_machine):

    print "exampletree_directive"

    env = state.document.settings.env
    suffix = env.config.source_suffix
    dirname = posixpath.dirname(env.docname)
    glob = 'glob' in options

    ret = []
    subnode = addnodes.toctree()
    includefiles = []
    includetitles = {}
    all_docnames = env.found_docs.copy()
    # don't add the currently visited file in catch-all patterns
    all_docnames.remove(env.docname)
    for entry in content:
        if not entry:
            continue

        if not glob:
            # look for explicit titles and documents ("Some Title <document>").
            m = caption_ref_re.match(entry)
            if m:
                docname = m.group(2)
                includetitles[docname] = m.group(1)
            else:
                docname = entry
            # remove suffixes (backwards compatibility)
#             if docname.endswith(suffix):
#                 docname = docname[:-len(suffix)]
            # absolutize filenames
            docname = posixpath.normpath(posixpath.join(dirname, docname))
#             if docname not in env.found_docs:
#                 ret.append(state.document.reporter.warning(
#                     'toctree references unknown document %r' % docname, line=lineno))
#             else:
            includefiles.append(docname)
        else:
            patname = posixpath.normpath(posixpath.join(dirname, entry))
            docnames = sorted(patfilter(all_docnames, patname))
            for docname in docnames:
                all_docnames.remove(docname) # don't include it again
                includefiles.append(docname)
            if not docnames:
                ret.append(state.document.reporter.warning(
                    'toctree glob pattern %r didn\'t match any documents' % entry,
                    line=lineno))
    subnode['includefiles'] = includefiles

    print "exampletree:", subnode

    subnode['includetitles'] = includetitles
    subnode['maxdepth'] = options.get('maxdepth', -1)
    subnode['glob'] = glob
    ret.append(subnode)

    print "exampletree:", subnode

    return ret

class ExampleTree(Directive):
    """
    Directive to notify Sphinx about the hierarchical structure of the docs,
    and to include a table-of-contents like tree in the current document.
    """

    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {
        'maxdepth': int,
        'glob': directives.flag,
        'hidden': directives.flag,
        'numbered': directives.flag,
    }

    def run(self):
        env = self.state.document.settings.env
        suffix = env.config.source_suffix
        glob = 'glob' in self.options

        ret = []
        # (title, ref) pairs, where ref may be a document, or an external link,
        # and title may be None if the document's title is to be used
        entries = []
        includefiles = []
        includetitles = {}
        all_docnames = env.found_docs.copy()
        # don't add the currently visited file in catch-all patterns
        all_docnames.remove(env.docname)
        for entry in self.content:
            if not entry:
                continue
            if not glob:
                # look for explicit titles ("Some Title <document>")
                m = caption_ref_re.match(entry)
                if m:
                    ref = m.group(2)
                    title = m.group(1)
                    docname = ref
                else:
                    ref = docname = entry
                    title = None
                # remove suffixes (backwards compatibility)
                if docname.endswith(suffix):
                    docname = docname[:-len(suffix)]
                # absolutize filenames
                docname = docname_join(env.docname, docname)
                if url_re.match(ref) or ref == 'self':
                    entries.append((title, ref))
#                 elif docname not in env.found_docs:
#                     ret.append(self.state.document.reporter.warning(
#                         'toctree references unknown document %r' % docname,
#                         line=self.lineno))
                else:
                    entries.append((title, docname))
                    includefiles.append(docname)
            else:
                patname = docname_join(env.docname, entry)
                docnames = sorted(patfilter(all_docnames, patname))
                for docname in docnames:
                    all_docnames.remove(docname) # don't include it again
                    entries.append((None, docname))
                    includefiles.append(docname)
                if not docnames:
                    ret.append(self.state.document.reporter.warning(
                        'toctree glob pattern %r didn\'t match any documents'
                        % entry, line=self.lineno))
        subnode = addnodes.toctree()
        subnode['parent'] = env.docname
        # entries contains all entries (self references, external links etc.)
        subnode['entries'] = entries
        # includefiles only entries that are documents
        subnode['includefiles'] = includefiles
        subnode['maxdepth'] = self.options.get('maxdepth', -1)
        subnode['glob'] = glob
        subnode['hidden'] = 'hidden' in self.options
        subnode['numbered'] = 'numbered' in self.options
        ret.append(subnode)
        return ret

def process_missing_example(app, env, node, contnode):
    print "app:", app
    print "env:", env
    print "node:", node
    print "contnode:", contnode

    return None

def doctree_read(app, doctree):
    print "doctree:", doctree.attributes['source']

def autodoc_process_docstring(app, what, name, obj, options, lines):
    if name.startswith("examples."):
        title = "Module :mod:`%s`" % name
        lines.insert(0, title)
        lines.insert(1, "-" * len(title))

def source_read(app, docname, source):
    print "source-read:", docname

def setup(app):
    app.add_directive('fipyexample', FiPyExample)
    app.add_directive('exampletree', ExampleTree)
#     app.add_directive('exampletree', exampletree_directive,
#                       content=1,
#                       arguments=(0, 0, 0),
#                       maxdepth=int,
#                       glob=directives.flag)
    app.connect('missing-reference', process_missing_example)
    app.connect('autodoc-process-docstring', autodoc_process_docstring)
    app.connect('doctree-read', doctree_read)
    app.connect('source-read', source_read)
