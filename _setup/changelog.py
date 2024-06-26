"""Uses basic authentication (Github username + password) to retrieve issues
from a repository that username has access to. Supports Github API v3.
Adapted from: https://gist.github.com/patrickfuller/e2ea8a94badc5b6967ef3ca0a9452a43
"""
from __future__ import print_function

import os
import textwrap
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
         "Only issues closed at or after this tag, SHA, or date are returned."),
        ('before=', None,
         "Only issues closed at or before this tag, SHA, or date are returned."),
        ('milestone=', None,
         "A string referring to a milestone by its title field. "
         "If the string `*` is passed, issues with any milestone are accepted. "
         "If the string `none` is passed, "
         "issues without milestones are returned. ")
     ]

    def initialize_options(self):
        import github

        self.repository = "usnistgov/fipy"
        self.tokenvar = "FIPY_GITHUB_TOKEN"
        self.username = None
        self.auth = None
        self.state = "closed"
        self.after = None
        self.before = None
        self.milestone = None

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

    def _printReST(self, issues, label):
        """Print section of issues to stdout
        """
        print()
        print(label)
        print("-" * len(label))
        print()

        for i, issue in issues.iterrows():
            print(issue.ReST)

    def _getMilestone(self, milestone):
        """Return Milestone with title of `milestone`

        If `milestone` is "*" or "none", returns unchanged
        """
        if milestone not in (None, "*", "none"):
            milestones = self.repo.get_milestones()
            milestones = [ms for ms in milestones if ms.title == self.milestone]
            try:
                milestone = milestones[0]
            except IndexError:
                raise KeyError("Milestone `{}` not found".format(self.milestone))

        return milestone

    def _getDateFromTagOrSHA(self, tagOrSHA):
        """Return date of tag or commit

        If tagOrSHA is None, return intact to allow all dates.
        If tagOrSHA corresponds to a tag.name, return date of its commit.
        If tagOrSHA corresponds to a commit SHA, return date of its commit.
        Else assume it's a date string.
        """
        if tagOrSHA is None:
            date = tagOrSHA
        else:
            tags = self.repo.get_tags()
            tags = [tag for tag in tags if tag.name == tagOrSHA]
            try:
                date = tags[0].commit.commit.author.date
            except IndexError:
                try:
                    commit = self.repo.get_commit(tagOrSHA)
                    date = commit.commit.author.date
                except:
                    date = tagOrSHA

        return date

    def read_pull(self, x):
        if x['pull_request'] is None:
            result = (False, None)
        else:
            pull = self.repo.get_pull(x.number)
            result = (pull.merged, pull.merged_at)
        return result

    def format_pull(self, x):
        prefix = "- "
        hang = " " * len(prefix)
        s = [textwrap.fill(x.title,
                           initial_indent=prefix,
                           subsequent_indent=hang)]
        s += [u"{}(`#{} <{}>`_)".format(hang,
                                        x.number,
                                        x.html_url)]
        s += ([u"{}Thanks to `@{} <{}>`_.".format(hang,
                                                  x.user.login,
                                                  x.user.html_url)]
               if x.user.login not in self.collaborators
               else [])
        return "\n".join(s)

    def format_issue(self, x):
        prefix = "- "
        hang = " " * len(prefix)
        s = [prefix + u"`#{} <{}>`_:".format(x.number,
                                             x.html_url)]
        s += [textwrap.fill(x.title,
                            initial_indent=hang,
                            subsequent_indent=hang)]
        return "\n".join(s)

    def run(self):
        """Requests issues from GitHub API and prints as ReST to stdout
        """
        import github
        import pandas as pd
        
        self.gh = github.Github(*self.auth)
        self.repo = self.gh.get_repo(self.repository)

        self.after = self._getDateFromTagOrSHA(self.after)
        self.before = self._getDateFromTagOrSHA(self.before)

        if self.after is not None:
            since = pd.to_datetime(self.after).to_pydatetime()
        else:
            since = github.GithubObject.NotSet

        milestone = self._getMilestone(self.milestone)
        if milestone is None:
            milestone = github.GithubObject.NotSet

        issues = self.repo.get_issues(state=self.state,
                                      since=since,
                                      milestone=milestone)
        self.collaborators = [c.login for c in self.repo.get_collaborators()]

        issues = [{
              'number': issue.number,
              'state': issue.state,
              'title': issue.title,
              'body': issue.body,
              'created_at': issue.created_at,
              'updated_at': issue.updated_at,
              'closed_at': issue.closed_at,
              'html_url': issue.html_url,
              'pull_request': issue.pull_request,
              'user': issue.user,
              'labels': [label.name for label in issue.labels]
            } for issue in issues]

        issues = pd.DataFrame(issues)

        issues = issues.sort_values(by=["number"], ascending=[False])
        wontfix = issues.labels.apply(lambda x: 'wontfix' in x)
        invalid = issues.labels.apply(lambda x: 'invalid' in x)
        question = issues.labels.apply(lambda x: 'question' in x)
        worksforme = issues.labels.apply(lambda x: 'worksforme' in x)
        duplicate = issues.labels.apply(lambda x: 'duplicate' in x)
        issues = issues[~wontfix & ~invalid & ~question & ~worksforme & ~duplicate]

        # fix the dates to reflect dates from original Trac issue tracker
        trac = (r" _Imported from trac ticket .*,  "
                r"created by .* on (.*), "
                r"last modified: (.*)_")
        olddates = issues.body.str.extract(trac).apply(pd.to_datetime)
        issues.loc[olddates[1].notna(), "created_at"] = olddates[0]
        issues.loc[olddates[1].notna(), "updated_at"] = olddates[1]
        issues.loc[((issues.state == "closed")
                    & olddates[1].notna()), "closed_at"] = olddates[1]

        if self.after is not None:
            issues = issues[issues['closed_at'] >= self.after]
        if self.before is not None:
            issues = issues[issues['closed_at'] <= self.before]

        ispull = issues['pull_request'].notna()
        isissue = ~ispull

        issues['merged'], issues['merged_at'] = zip(*issues.apply(self.read_pull, axis=1))

        issues.loc[ispull, 'ReST'] = issues.apply(self.format_pull, axis=1)

        issues.loc[isissue, 'ReST'] = issues.apply(self.format_issue, axis=1)

        self._printReST(issues[ispull & issues['merged']], "Pulls")
        self._printReST(issues[isissue], "Fixes")
