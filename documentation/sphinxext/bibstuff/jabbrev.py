#! /usr/bin/env python
# File: jabbrev.py

"""
A utility to read journal abbreviations from a text file and modify a
bibtex database to use the abbreviations.  Right now it just replaces
long versions with abbreviations in bib file.

:author: Dylan Schwilk
:contact: http://www.pricklysoft.org 
:copyright: 2006 by Dylan Schwilk and Alan G Isaac
:license: MIT (see `license.txt`_)
:date: 2006-09-15
:TODO: Move this to a scripts directory or the examples directory?

.. _license.txt: ./license.txt
"""

__docformat__ = "restructuredtext en"
__version__ = "1.5"
__author__ = "Dylan W. Schwilk"

import string, sys, re
from simpleparse.parser import Parser
from simpleparse.dispatchprocessor import * 

import bibfile,  bibgrammar


def MakeMap(jlist):
    """makes a dictionary of journals.
    key long name by default, data abbreviation.
    """
    jmap = {}
    for j in jlist:
        l = j.split("=")
        jmap[l[1].strip()] = l[0].strip()
    return jmap

def Translate(bib, jmap):
    """Translate"""
    for entry in bib.entries:
        j = entry['journal']
        if j :
            entry['journal'] = jmap.get(j, j)
        print entry   
            
      
def main():
    '''Command-line tool.  See jabbrev.py -h for help'''

    input = sys.stdin
    output = sys.stdout
    
    try:
        from optparse import OptionParser
    except (ImportError, AttributeError):
        try:
            from optik import OptionParser
        except (ImportError, AttributeError):
            print "jabbrev needs python 2.3 or Greg Ward's optik module."
    
    usage = "usage: %prog [options] DATABASE_FILE [ABBREVIATION_FILE] "
    parser = OptionParser(usage=usage, version ="%prog " + __version__)
     
    (options, args) = parser.parse_args()

    if len(args) > 2 :
        print "Too many arguments"
        sys.exit(1)
    try :
        # bibtex file is last argument
        bib = open(args[0]).read()
    except :
        print "No bibtex file found."
        sys.exit(1)
    if len(args) == 2 :
        input = open(args[1])
    else :
        input = sys.stdin.read()
              
    bfile = bibfile.BibFile()
    bibgrammar.Parse(bib, bfile)

    journals = MakeMap(input.readlines())    
    Translate(bfile, journals)

if __name__ == '__main__':
    main()
