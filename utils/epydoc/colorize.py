#
# epydoc.html: HTML colorizers
# Edward Loper
#
# Created [10/16/02 09:49 PM]
# $Id$
#

r"""
Functions to produce colorized HTML code for various objects.
Currently, C{colorize} defines functions to colorize regular
expressions and doctest blocks.

@var RE_TAG: The CSS class for colorizing regular expressions.

@var ANY_TAG: The CSS class for colorizing C{"."} in regular
expressions.

@var ESCAPE_TAG: The CSS class for colorizing escaped characters (such
as C{r"\("}) in regular expressions.
    
@var CATEGORY_TAG: The CSS class for colorizing character categories
(such as C{r"\d"})) in regular expressions.
    
@var AT_TAG: The CSS class for colorizing character locations (such as
C{"^"}) in regular expressions.
    
@var BRANCH_TAG: The CSS class for colorizing C{"|"} in regular
expressions.
    
@var STAR_TAG: The CSS class for colorizing C{"*"} and C{"*?"} in
regular expressions.
    
@var PLUS_TAG: The CSS class for colorizing C{"+"} and C{"+?"} in
regular expressions.
    
@var QMRK_TAG: The CSS class for colorizing C{"?"} and C{"??"} in
regular expressions.
    
@var RNG_TAG: The CSS class for colorizing repeat ranges (such as
C{"a{3,8}"}) in regular expressions.
    
@var PAREN_TAG: The CSS class for colorizing parenthases in regular
expressions.
    
@var CHOICE_TAG: The CSS class for colorizing character choice
expressions (such as C{"[abc]"}) in regular expressions.
    
@var ASSERT_TAG: The CSS class for colorizing assertions (such as
C{"(?=abc)"}) in regular expressions.
    
@var REF_TAG: The CSS class for colorizing references (such as
C{r"\1"}) in regular expressions.

@var _PROMPT_RE: The regular expression used to find Python prompts
(">>>" and "...") in doctest blocks.

@var _DOCTEST_RE: The regular expression used by L{_doctest_sub} to
colorize doctest blocks.
"""
__docformat__ = 'epytext en'

import sys, sre_parse, sre, re
import sre_constants

##################################################
## Regular expression colorizer
##################################################

# HTML tags for colorize_re
RE_TAG         = 're'
ANY_TAG        = 're-char'      # r"."
ESCAPE_TAG     = 're-char'      # r"\("
CATEGORY_TAG   = 're-char'      # r"\d"
AT_TAG         = 're-char'      # r"^"
BRANCH_TAG     = 're-op'        # r"a|b|c"
STAR_TAG       = 're-op'        # r"a*"
PLUS_TAG       = 're-op'        # r"a+"
QMRK_TAG       = 're-op'        # r"a?"
RNG_TAG        = 're-op'        # r"a{3,8}"
PAREN_TAG      = 're-group'     # r"(abc)"
CHOICE_TAG     = 're-group'     # r"[abc]"
ASSERT_TAG     = 're-group'     # r"(?=foo)"
REF_TAG        = 're-ref'       # r"\1"

def colorize_re(regexp):
    r"""
    @return: The HTML code for a colorized version of the pattern for
        the given SRE regular expression.  If C{colorize_re} can't
        figure out how to colorize the regexp, then it will simply return
        the (uncolorized) pattern, with C{'&'}, C{'<'}, and C{'>'}
        escaped as HTML entities.  The colorized expression includes
        spans with the following css classes:
          - X{re}: The entire regular expression.
          - X{re-char}: Special characters (such as C{'.'}, C{'\('}), 
            character categories (such as C{'\w'}), and locations
            (such as C{'\b'}).
          - X{re-op}: Operators (such as C{'*'} and C{'|'}).
          - X{re-group}: Grouping constructs (such as C{'(...)'}).
          - X{re-ref} References (such as C{'\1'})
    @rtype: C{string}
    @param regexp: The regular expression to colorize.
    @type regexp: C{SRE_Pattern} or C{string}
    """
    try:
        if type(regexp) == type(''): regexp = sre.compile(regexp)
        tree = sre_parse.parse(regexp.pattern, regexp.flags)
        return ('<span class="%s">%s</span>' %
                (RE_TAG, _colorize_re(tree, 1)))
    except:
        try:
            pat = regexp.pattern
            pat = pat.replace('&', '&amp;')
            pat = pat.replace('<', '&lt;')
            pat = pat.replace('>', '&gt;')
            return '<span class="%s">%s</span>' % (RE_TAG, pat)
        except:
            try:
                str = `regexp`
                str = str.replace('&', '&amp;')
                str = str.replace('<', '&lt;')
                str = str.replace('>', '&gt;')
                return str
            except: return '<span class="%s">...</span>' % RE_TAG
    
def _colorize_re(tree, noparen=0):
    """
    Recursively descend the given regexp parse tree to produce the
    HTML code for a colorized version of the regexp.

    @param tree: The regexp parse tree for the regexp that should be
        colorized.
    @type tree: L{sre_parse.SubPattern}
    @param noparen: If true, then don't include parenthases around the
        expression in C{tree}, even if it contains multiple elements.
    @type noparen: C{boolean}
    @return: The HTML code for a colorized version of C{tree}
    @rtype: C{string}
    """
    str = ''
    if len(tree) > 1 and not noparen:
        str += '<span class="%s">(</span>' % PAREN_TAG
    for elt in tree:
        op = elt[0]
        args = elt[1]

        if op == sre_constants.LITERAL:
            c = chr(args)
            if c == '&': str += '&amp;'
            elif c == '<': str += '&lt;'
            elif c == '>': str += '&gt;'
            elif c == '\t': str += r'<span class="%s">\t</span>' % ESCAPE_TAG
            elif c == '\n': str += r'<span class="%s">\n</span>' % ESCAPE_TAG
            elif c == '\r': str += r'<span class="%s">\r</span>' % ESCAPE_TAG
            elif c == '\f': str += r'<span class="%s">\f</span>' % ESCAPE_TAG
            elif c == '\v': str += r'<span class="%s">\v</span>' % ESCAPE_TAG
            elif c in '.^$\\*+?{}[]|()':
                str += '<span class="%s">\\%c</span>' % (ESCAPE_TAG, c)
            else: str += chr(args)
            continue
        
        elif op == sre_constants.ANY:
            str += '<span class="%s">.</span>' % ANY_TAG
            
        elif op == sre_constants.BRANCH:
            if args[0] is not None:
                raise ValueError('Branch expected None arg but got %s'
                                 % args[0])
            VBAR = '<span class="%s">|</span>' % BRANCH_TAG
            str += VBAR.join([_colorize_re(item,1) for item in args[1]])
            
        elif op == sre_constants.IN:
            if (len(args) == 1 and args[0][0] == sre_constants.CATEGORY):
                str += _colorize_re(args)
            else:
                str += '<span class="%s">[</span>' % CHOICE_TAG
                str += _colorize_re(args, 1)
                str += '<span class="%s">]</span>' % CHOICE_TAG
                
        elif op == sre_constants.CATEGORY:
            str += '<span class="%s">' % CATEGORY_TAG
            if args == sre_constants.CATEGORY_DIGIT: str += r'\d'
            elif args == sre_constants.CATEGORY_NOT_DIGIT: str += r'\D'
            elif args == sre_constants.CATEGORY_SPACE: str += r'\s'
            elif args == sre_constants.CATEGORY_NOT_SPACE: str += r'\S'
            elif args == sre_constants.CATEGORY_WORD: str += r'\w'
            elif args == sre_constants.CATEGORY_NOT_WORD: str += r'\W'
            else: raise ValueError('Unknown category %s' % args)
            str += '</span>'
            
        elif op == sre_constants.AT:
            str += '<span class="%s">' % AT_TAG
            if args == sre_constants.AT_BEGINNING_STRING: str += r'\A'
            elif args == sre_constants.AT_BEGINNING: str += r'^'
            elif args == sre_constants.AT_END: str += r'$'
            elif args == sre_constants.AT_BOUNDARY: str += r'\b'
            elif args == sre_constants.AT_NON_BOUNDARY: str += r'\B'
            elif args == sre_constants.AT_END_STRING: str += r'\Z'
            else: raise ValueError('Unknown position %s' % args)
            str += '</span>'
            
        elif op == sre_constants.MAX_REPEAT:
            min = args[0]
            max = args[1]
            if max == sre_constants.MAXREPEAT:
                if min == 0:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">*</span>' % STAR_TAG
                elif min == 1:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">+</span>' % PLUS_TAG
                else:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">{%d,}</span>' % (RNG_TAG, min)
            elif min == 0:
                if max == 1:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">?</span>' % QMRK_TAG
                else:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">{,%d}</span>' % (RNG_TAG, max)
            elif min == max:
                str += _colorize_re(args[2])
                str += '<span class="%s">{%d}</span>' % (RNG_TAG, max)
            else:
                str += _colorize_re(args[2])
                str += '<span class="%s">{%d,%d}</span>' % (RNG_TAG, min, max)

        elif op == sre_constants.MIN_REPEAT:
            min = args[0]
            max = args[1]
            if max == sre_constants.MAXREPEAT:
                if min == 0:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">*?</span>' % STAR_TAG
                elif min == 1:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">+?</span>' % PLUS_TAG
                else:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">{%d,}?</span>' % (RNG_TAG, min)
            elif min == 0:
                if max == 1:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">??</span>' % QMRK_TAG
                else:
                    str += _colorize_re(args[2])
                    str += '<span class="%s">{,%d}?</span>' % (RNG_TAG, max)
            elif min == max:
                str += _colorize_re(args[2])
                str += '<span class="%s">{%d}?</span>' % (RNG_TAG, max)
            else:
                str += _colorize_re(args[2])
                str += '<span class="%s">{%d,%d}?</span>'%(RNG_TAG, min, max)

        elif op == sre_constants.SUBPATTERN:
            if args[0] is None:
                str += '<span class="%s">(?:</span>' % PAREN_TAG
            elif type(args[0]) == type(0):
                # This is cheating:
                str += '<span class="%s">(</span>' % PAREN_TAG
            else:
                str += '<span class="%s">(?P&lt;</span>' % PAREN_TAG
                str += '<span class="%s">%s</span>' % (REF_TAG, args[0])
                str += '<span class="%s">&gt;</span>' % PAREN_TAG
            str += _colorize_re(args[1], 1)
            str += '<span class="%s">)</span>' % PAREN_TAG

        elif op == sre_constants.GROUPREF:
            str += '<span class="%s">\\%d</span>' % (REF_TAG, args)

        elif op == sre_constants.RANGE:
            str += ('%c<span class="%s">-</span>%c' %
                    (chr(args[0]), CHOICE_TAG, chr(args[1])))
            
        elif op == sre_constants.NEGATE:
            str += '<span class="%s">^</span>' % CHOICE_TAG

        elif op == sre_constants.ASSERT:
            if args[0]: str += '<span class="%s">(?=</span>' % ASSERT_TAG
            else: str += '<span class="%s">(?&lt;=</span>' % ASSERT_TAG
            str += ''.join(_colorize_re(args[1], 1))
            str += '<span class="%s">)</span>' % ASSERT_TAG
                           
        elif op == sre_constants.ASSERT_NOT:
            if args[0]: str += '<span class="%s">(?!</span>' % ASSERT_TAG
            else: str += '<span class="%s">(?&lt;!</span>' % ASSERT_TAG
            str += ''.join(_colorize_re(args[1], 1))
            str += '<span class="%s">)</span>' % ASSERT_TAG

        elif op == sre_constants.NOT_LITERAL:
            lit = _colorize_re( ((sre_constants.LITERAL, args),) )
            str += ('<span class="%s">[^</span>%s<span class="%s">]</span>' %
                    (CHOICE_TAG, lit, CHOICE_TAG))
        else:
            print 'UNKNOWN ELT', elt[0], elt
    if len(tree) > 1 and not noparen: 
        str += '<span class="%s">)</span>' % PAREN_TAG
    return str

##################################################
## Doctest block colorizer
##################################################

# Regular expressions for colorize_doctestblock
_KEYWORDS = ["del", "from", "lambda", "return", "and", "or", "is", 
             "global", "not", "try", "break", "else", "if", "elif", 
             "while", "class", "except", "import", "pass", "raise",
             "continue", "finally", "in", "print", "def", "for"]
_KEYWORD = '|'.join([r'(\b%s\b)' % _KW for _KW in _KEYWORDS])
_STRING = '|'.join([r'("""("""|.*?((?!").)"""))', r'("("|.*?((?!").)"))',
                    r"('''('''|.*?[^\\']'''))", r"('('|.*?[^\\']'))"])
_STRING = _STRING.replace('"', '&quot;') # Careful with this!
_COMMENT = '(#.*?$)'
_PROMPT = r'(^\s*(&gt;&gt;&gt;|\.\.\.)(\s|$))'
_PROMPT_RE = re.compile(_PROMPT, re.MULTILINE | re.DOTALL)
_DOCTEST_RE = re.compile('|'.join([_STRING, _COMMENT, _KEYWORD]),
                          re.MULTILINE | re.DOTALL)
del _KEYWORDS, _KEYWORD, _STRING, _COMMENT, _PROMPT, _KW

def colorize_doctestblock(str):
    """
    @return: The HTML code for a colorized version of a given doctest
        block.  In particular, this identifies spans with the
        following css classes:
          - X{py-src}: The Python source code.
          - X{py-prompt}: The ">>>" and "..." prompts.
          - X{py-string}: Strings in the Python source code.
          - X{py-comment}: Comments in the Python source code.
          - X{py-keyword}: Keywords in the Python source code.
          - X{py-output}: Python's output (lines without a prompt).
        The string that is passed to colorize_doctest should already
        have HTML characters escaped (e.g., C{">"} should be encoded
        as C{"&gt;"}).
    @type str: C{string}
    @param str: The contents of the doctest block to be colorized.
    @rtype: C{string}
    """
    pysrc = pyout = ''
    outstr = ''
    for line in str.split('\n')+['\n']:
        if _PROMPT_RE.match(line):
            if pyout:
                outstr += ('<span class="py-output">%s</span>\n\n' %
                           pyout.strip())
                pyout = ''
            pysrc += line+'\n'
        else:
            if pysrc:
                # Prompt over-rides other colors (incl string)
                pysrc = _DOCTEST_RE.sub(_doctest_sub, pysrc)
                pysrc = _PROMPT_RE.sub(r'<span class="py-prompt">'+
                                       r'\1</span>', pysrc)
                outstr += ('<span class="py-src">%s</span>\n'
                           % pysrc.strip())
                pysrc = ''
            pyout += line+'\n'
    if pyout.strip():
        outstr += ('<span class="py-output">%s</span>\n' %
                   pyout.strip())
    return outstr.strip()
    
def _doctest_sub(match):
    """
    This helper function is used by L{colorize_doctestblock} to
    add colorization to matching expressions.  It is called by
    C{_DOCTEST_RE.sub} with an expression that matches
    C{_DOCTEST_RE}.

    @return: The HTML code for the colorized expression.
    @rtype: C{string}
    @see: L{_DOCTEST_RE}
    """
    str = match.group()
    if str[:1] == "'" or str[:6] == '&quot;':
        return '<span class="py-string">%s</span>' % str
    elif str[:1] in '#':
        return '<span class="py-comment">%s</span>' % str
    else:
        return '<span class="py-keyword">%s</span>' % str

