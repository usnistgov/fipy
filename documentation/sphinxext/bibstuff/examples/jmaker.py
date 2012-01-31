#!usr/bin/python

# simple script to take journal list from
# http://www.csa.com/htbin/sjldisp.cgi?filename=/wais/data/srcjnl/biologset
# and produce a format readable by jabbrev.py


import sys, string

input = sys.stdin.readlines()

for line in input:
    names = line.split(', (')
    journal = names[0].strip()
    try :
        info = names[1].split(')')
        abbrev = info[0].strip()
    except :
        abbrev = journal
    print "%s = %s" % (abbrev, journal)
