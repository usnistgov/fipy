## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "Change_log.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

"""Uses basic authentication (Github username + password) to retrieve issues
from a repository that username has access to. Supports Github API v3.
Forked from: patrickfuller/github_issues_to_csv.py
"""

from distutils.core import Command

__all__ = ["Change_log"]

class Change_log(Command):
    description = "Generate ReST change log from github issues and pull requests"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [
        ('repository=', None,
         "GitHub repository to obtain issues from (default: 'usnistgov/fipy')"),
        ('username=', None,
         "GitHub username to authenticate as (default: None). Note: GitHub limits the rate of unauthenticated queries: https://developer.github.com/v3/#rate-limiting")
        ('state=', None,
         "Indicates the state of the issues to return. Can be either `open`, `closed`, or `all`. (default: `closed`)"),
     ]

    def initialize_options(self):
        self.repository = "usnistgov/fipy"
        self.username = None
        self.auth = None
        self.state = "closed"

    def finalize_options(self):
        if self.username is not None:
            from getpass import getpass
            
            password = getpass("Password for 'https://{}@github.com': ".format(self.username))
            self.auth = (username, password)

    def run(self):
        """Requests issues from GitHub API and writes to JSON file.
        """
        url = 'https://api.github.com/repos/{}/issues?state={}'.format(name, state)
        r = requests.get(url, auth=auth)
        
        jsondat = []
        
        while True:
            if r.status_code != 200:
                raise Exception(r.status_code)

            jsondat += r.json()
            
            if 'link' not in r.headers:
                break
            
            pages = {rel[6:-1]: url[url.index('<')+1:-1] for url, rel in
                     (link.split(';') for link in
                      r.headers['link'].split(','))}

            if 'next' not in pages:
                break
                
            r = requests.get(pages['next'], auth=auth)
            
        jsonfilename = '{}-issues.json'.format(name.replace('/', '-'))
            
        with open(jsonfilename, 'w') as jsonfile:
            json.dump(jsondat, jsonfile, indent=2)
