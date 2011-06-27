##None of this is relevant as it's been merged directly into efficienct_test.py

import urllib, urllib2
              
CODESPEED_URL = "http://localhost:8000/"
                   
def add(data):
    print data
    params = urllib.urlencode(data)
    response = "None"
    print "Executable %s, revision %s, benchmark %s" % (data['executable'],data['commitid'], data['benchmark'])
    f = urllib2.urlopen(CODESPEED_URL + 'result/add/', params)  
#    print type(f) , '\n' , f
    response = f.read()
    f.close()
    print "Server (%s) response: %s\n" % (CODESPEED_URL, response)
                       
if __name__ == '__main__':
    from datetime import datetime
    data = {
        'commitid': '6',
        'branch': 'default',#Always use default for trunk/master/tip
        'project': 'Prototype Test',
        'revision_date': '', # Optional. Default is taken either
        # from VCS integration or from current date
        'executable': 'datatest.py',
        'benchmark': 'float',
        'environment': "Test",
        'result_value': 0.14,
        'result_date': datetime.today(), # Optional
        #    'std_dev': 1.11111, # Optional. Default is blank
        #    'max': 4001.6, # Optional. Default is blank
        #    'min': 3995.1, # Optional. Default is blank
        }

    add(data)
