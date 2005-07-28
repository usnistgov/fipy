#
# test/epytext.py: tests for epytext markup language
# Edward Loper
#
# Created [10/30/02 12:06 PM]
# $Id$
#

"""
Regression testing for the epytext markup language.

The test cases currently implemented are by no means comprehensive.
"""
__docformat__ = 'epytext en'

if __name__ == '__main__':
    import epydoc.markup.epytext; reload(epydoc.markup.epytext)

import unittest, re
from epydoc.markup.epytext import *

##//////////////////////////////////////////////////////
##  Parse Tests
##//////////////////////////////////////////////////////

class ParseTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def failIfParseError(self, text, errors):
        estr = '\n' + '~'*70 + '\n' + text
        if errors:
            estr += '\n'+'_'*70+'\nERRORS:\n'
            for e in errors: estr += '%s\n' % e
        if errors:
            self.fail(estr+'~'*70)

    def checkParse(self, epytext, xml=None):
        """
        Parse C{epytext}, and check that it generates xml output
        X{xml}, with no warnings or errors.
        """
        errors = []
        out = parse(epytext, errors)
        if out is None: out = ''
        else: out = out.childNodes[0].toxml().strip()
        if out[:9] == '<epytext>' and out[-10:] == '</epytext>':
            out = out[9:-10]
            
        self.failIfParseError(epytext, errors)
        if xml:
            self.failUnlessEqual(`out`, `xml.strip()`)

    def checkParseError(self, epytext, errtype, linenum):
        errors = []
        out = parse(epytext, errors)

        for err in errors:
            if isinstance(err, errtype) and err.linenum() == linenum:
                errors.remove(err)
                break
        else:
            self.fail("No %s generated on line %s" %
                      (errtype.__name__, linenum))

        self.failIfParseError(epytext, errors)

    def testPara(self):
        self.checkParse("""
        this is one paragraph.

        This is
        another.

        This is a third.""",
                        '<para>this is one paragraph.</para>'+
                        '<para>This is another.</para><para>This is a '+
                        'third.</para>')

    def testUnindentedFields(self):
        """
        Make sure that unindented fields are allowed.
        """
        self.checkParse("""
        This is a paragraph.
        
        @foo: This is a field.""")
        
        self.checkParse("""
        This is a paragraph.
        @foo: This is a field.""")
        
        self.checkParse("""
        This is a paragraph.
           @foo: This is a field.
             Hello.""")
        
        self.checkParse("""
        This is a paragraph.
           @foo: This is a field.
             Hello.""")
        self.checkParse("""Paragraph\n@foo: field""")
        self.checkParse("""Paragraph\n\n@foo: field""")
        self.checkParse("""\nParagraph\n@foo: field""")

    def testUnindentedList(self):
        """
        Make sure that unindented lists are not allowed.
        """
        self.checkParseError("""
        This is a paragraph.
        
        - This is a list item.""", StructuringError, 4)
        
        self.checkParseError("""
        This is a paragraph.
        - This is a list item.""", StructuringError, 3)
        
        self.checkParseError("""
        This is a paragraph.
           - This is a list item.
             Hello.
             - Sublist item.""", StructuringError, 5)
        
        self.checkParseError("""
        This is a paragraph.
           - This is a list item.
             Hello.
             
             - Sublist item.""", StructuringError, 6)
        self.checkParseError("""\nParagraph\n- list item""",
                             StructuringError, 3)
        self.checkParseError("""Paragraph\n\n- list item""",
                             StructuringError, 3)
        self.checkParseError("""\nParagraph\n- list item""",
                             StructuringError, 3)

        # Special case if there's text on the same line as the opening
        # quote..
        self.checkParse("""Paragraph\n- list item""",
                        "<para>Paragraph</para><ulist>"+
                        "<li><para>list item</para></li></ulist>")

    def testIndentedList(self):
        """
        Make sure that indented lists are allowed.
        """
        list1 = ("<para>This is a paragraph.</para><ulist>"+
                 "<li><para>This is a list item.</para></li>"+
                 "</ulist><para>This is a paragraph</para>")
        list2 = '<ulist><li><para>This is a list item.</para></li></ulist>'
        
        self.checkParse('This is a paragraph.\n  - This is a list item.\n'+
                        'This is a paragraph', list1)
        self.checkParse('This is a paragraph.\n\n  - This is a list item.'+
                        '\n\nThis is a paragraph', list1)
        self.checkParse("""
        This is a paragraph.
        
          - This is a list item.
          
        This is a paragraph""", list1)
        self.checkParse("""
        This is a paragraph.
        
              - This is a list item.
        This is a paragraph""", list1)
        self.checkParse("""
        - This is a list item.""", list2)
        self.checkParse("""- This is a list item.""", list2)
        self.checkParse("""\n- This is a list item.""", list2)

    def testListBasic(self):
        P1 = "This is a paragraph."
        P2 = "This is a \nparagraph."
        LI1 = "  - This is a list item."
        LI2 = "\n  - This is a list item."
        LI3 = "  - This is a list\n  item."
        LI4 = "\n  - This is a list\n  item."
        PARA = ('<para>This is a paragraph.</para>')
        ONELIST = ('<ulist><li><para>This is a '+
                   'list item.</para></li></ulist>')
        TWOLIST = ('<ulist><li><para>This is a '+
                   'list item.</para></li><li><para>This is a '+
                   'list item.</para></li></ulist>')

        for p in (P1, P2):
            for li1 in (LI1, LI2, LI3, LI4):
                self.checkParse(li1, ONELIST)
                self.checkParse('%s\n%s' % (p, li1), PARA+ONELIST)
                self.checkParse('%s\n%s' % (li1, p), ONELIST+PARA)
                self.checkParse('%s\n%s\n%s' % (p, li1, p),
                                PARA+ONELIST+PARA)
            
                for li2 in (LI1, LI2, LI3, LI4):
                    self.checkParse('%s\n%s' % (li1, li2), TWOLIST)
                    self.checkParse('%s\n%s\n%s' % (p, li1, li2), PARA+TWOLIST)
                    self.checkParse('%s\n%s\n%s' % (li1, li2, p), TWOLIST+PARA)
                    self.checkParse('%s\n%s\n%s\n%s' % (p, li1, li2, p),
                                    PARA+TWOLIST+PARA)

        LI5 = "  - This is a list item.\n\n    It contains two paragraphs."
        LI5LIST = ('<ulist><li><para>This is a list item.</para>'+
                   '<para>It contains two paragraphs.</para></li></ulist>')
        self.checkParse(LI5, LI5LIST)
        self.checkParse('%s\n%s' % (P1, LI5), PARA+LI5LIST)
        self.checkParse('%s\n%s\n%s' % (P2, LI5, P1), PARA+LI5LIST+PARA)

        LI6 = ("  - This is a list item with a literal block::\n" +
               "    hello\n      there")
        LI6LIST = ('<ulist><li><para>This is a list item with a literal '+
                   'block:</para><literalblock> hello\n   there'+
                   '</literalblock></li></ulist>')
        self.checkParse(LI6, LI6LIST)
        self.checkParse('%s\n%s' % (P1, LI6), PARA+LI6LIST)
        self.checkParse('%s\n%s\n%s' % (P2, LI6, P1), PARA+LI6LIST+PARA)

    def testListItemWrap(self):
        LI = "- This is a list\n  item."
        ONELIST = ('<ulist><li><para>This is a '+
                   'list item.</para></li></ulist>')
        TWOLIST = ('<ulist><li><para>This is a '+
                   'list item.</para></li><li><para>This is a '+
                   'list item.</para></li></ulist>')
        for indent in ('', '  '):
            for nl1 in ('', '\n'):
                self.checkParse(nl1+indent+LI, ONELIST)
                for nl2 in ('\n', '\n\n'):
                    self.checkParse(nl1+indent+LI+nl2+indent+LI, TWOLIST)
            

##//////////////////////////////////////////////////////
##  Test Suite & Test Running
##//////////////////////////////////////////////////////

def testsuite():
    """
    Return a PyUnit testsuite for the epytext module.
    """
    
    tests = unittest.TestSuite()

    parse_tests = unittest.makeSuite(ParseTestCase, 'test')
    tests = unittest.TestSuite( (tests, parse_tests) )

    return tests

def test():
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(testsuite())

if __name__ == '__main__':
    test()

