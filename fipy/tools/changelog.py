## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "changelog.py"
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
Adapted from: https://gist.github.com/patrickfuller/e2ea8a94badc5b6967ef3ca0a9452a43
"""

import os
from distutils.core import Command

__all__ = ["changelog"]

class changelog(Command):
    description = "Generate ReST change log from github issues and pull requests"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [
        ('repository=', None,
         "GitHub repository to obtain issues from (default: 'usnistgov/fipy')"),
        ('tokenvar=', None,
         "Environment variable holding GitHub personal access token "
         "with 'repo' scope (default: 'FIPY_GITHUB_TOKEN')"),
        ('username=', None,
         "GitHub username to authenticate as (default: None). "
         "Supersedes `tokenvar`. "
         "Note: GitHub limits the rate of unauthenticated queries: "
         "https://developer.github.com/v3/#rate-limiting"),
        ('state=', None,
         "Indicates the state of the issues to return. "
         "Can be either `open`, `closed`, or `all`. (default: `closed`)"),
        ('after=', None,
         "Only issues closed at or after this tag are returned."),
        ('before=', None,
         "Only issues closed at or before this tag are returned.")
     ]

    def initialize_options(self):
        self.repository = "usnistgov/fipy"
        self.tokenvar = "FIPY_GITHUB_TOKEN"
        self.username = None
        self.auth = None
        self.state = "closed"
        self.after = None
        self.before = None

    def finalize_options(self):
        if self.username is not None:
            from getpass import getpass

            password = getpass("Password for 'https://{}@github.com': ".format(self.username))
            self.auth = (username, password)
        else:
            try:
                self.auth = (os.environ[self.tokenvar],)
            except KeyError:
                pass

    def run(self):
        """Requests issues from GitHub API and prints as ReST to stdout
        """
        import github
        import pandas as pd
        
        g = github.Github(*self.auth)
        repo = g.get_repo(self.repository)
        
        issues = [{
              'number': issue.number,
              'state': issue.state,
              'title': issue.title,
              'body': issue.body,
              'created_at': issue.created_at,
              'updated_at': issue.updated_at,
              'closed_at': issue.closed_at,
              'url': issue.url,
              'pull_request': issue.pull_request,
              'user': issue.user
            } for issue in repo.get_issues(state=self.state)]
            
        issues = pd.DataFrame(issues)

        if self.after is not None:
            issues = issues[issues['closed_at'] >= self.after]
        if self.before is not None:
            issues = issues[issues['closed_at'] <= self.before]

        pulls = issues[issues['pull_request'].notna()]
        issues = issues[issues['pull_request'].isna()]

        print issues
