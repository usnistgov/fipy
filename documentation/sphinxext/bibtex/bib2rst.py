#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
convert bibtex files to ReSt documents

It can also parse the occuring citations and include
only the relevant ones into the reference file.
"""
import subprocess
import os
import fnmatch
import sys

#import dependencies
import simpleparse

#local imports

from bibstuff import bibfile, bibgrammar, bibstyles
from bibstuff.bibstyles import default as style
from bibstuff import ebnf_sp
from bibstuff import bib4txt

def convert_all_bib(input_bib_directory, bibstuff_path):
    """Sphinx extension to convert BibTeX file to ReSt files.
    
    """
    #print  'CURR:', os.path.abspath(os.path.curdir)
    input_bib_directory = os.path.abspath(input_bib_directory)
    bibstuff_path = os.path.abspath(bibstuff_path)
    sys.path.append(bibstuff_path)
    for name in fnmatch.filter(os.listdir(input_bib_directory),'*.bib'):
#            print fnmatch.filter(os.listdir(input_directory),'*.bib')
            bib_in = os.path.join(input_bib_directory, name)
            # TODO: change into logging
#            print bib_in 
            #bib_in= bibfile_name = '_static/example.bib'
            bib_rst_out= bib_in+'_all.txt'
#            print bib_rst_out
#            print bibstuff_path
#            bibstuff_converter = os.path.join(bibstuff_path, 'bib4txt.py')
##            print bibstuff_converter
#            subprocess.call(['python', bibstuff_converter, '--all', '-no '+bib_rst_out,
#                                        bib_in])
            
            
    bibfile_processor = bibfile.BibFile()
    bibfile_as_string = open(bib_in,'r').read()
    bibgrammar.Parse(bibfile_as_string, bibfile_processor)
    entries = bibfile_processor.entries
    parsed_bibfile=bibfile_processor
    citation_manager = style.CitationManager([parsed_bibfile], keys=None,
                                             citation_template=style.CITATION_TEMPLATE)
    citation_manager.make_citations(entries)
    myres = citation_manager.make_citations(entries)
    out = open(bib_rst_out, 'w')
    out.write(myres)
    out.close()
    
def convert_only_cited(input_rst_directory, input_bib_path,
                                    source_suffix, bibstuff_path):
#    print  input_rst_directory
    input_rst_directory = os.path.abspath(input_rst_directory)
#     print  input_rst_directory
    
#    print input_bib_path
    input_bib_path = os.path.abspath(input_bib_path)
#    print input_bib_path
    bibstuff_path = os.path.abspath(bibstuff_path)
#    print  bibstuff_path
    sys.path.append(bibstuff_path)
    output_bib_file = input_bib_path+'_cited.txt'
#    print output_bib_file
    bibstuff_converter = os.path.join(bibstuff_path, 'bib4txt.py')
    src_as_string = ''
#    now add all strings from all files
    #for file in fnmatch.filter(os.walk(input_rst_directory),'*.rst'):
    ### http://www.brunningonline.net/simon/blog/archives/002022.html
    for path, dirs, files in os.walk(input_rst_directory): #os.walk(os.getcwd()):
        for file in [os.path.abspath(os.path.join(path, filename))
                         for filename in files if fnmatch.fnmatch(filename, '*'+source_suffix)]:
            #file = fnmatch.filter(file,'*.rst')
            print file
            rst_in = os.path.abspath(file)
            rst_reader = open(rst_in, 'r')
            rst_as_string = rst_reader.read()
            rst_reader.close()
            src_as_string = src_as_string+rst_as_string
    
#    print src_as_string
    
    ebnf_dec = ebnf_sp.cites_rest
    src_parser = cite_parser = simpleparse.parser.Parser(ebnf_dec, root='src')
    
    
    bibfile_processor = bibfile.BibFile()
    bibfile_as_string = open(input_bib_path,'r').read()
    bibgrammar.Parse(bibfile_as_string, bibfile_processor)
    parsed_bibfile = bibfile_processor
    
    result = bib4txt.make_text_output(src_as_string,
                     src_parser,
                     parsed_bibfile,
                     style,
                     citations_only=True)
    
    output = open(output_bib_file,'w')
    output.write(result)        
    output.close()
            
        
        