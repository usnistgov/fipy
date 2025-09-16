"""Uses personal access token or basic authentication (Github username +
password) to retrieve issues from a repository that username has access to.
Supports Github API version 3.
Adapted from: https://gist.github.com/patrickfuller/e2ea8a94badc5b6967ef3ca0a9452a43
"""

import textwrap
import typer
from typer import Argument, Option
from typing_extensions import Annotated

__all__ = ["app"]

class ChangeLog(object):
    """Generate reST change log from github issues and pull requests
    """

    def __init__(self,
        repository,
        auth,
        username,
        state,
        after,
        before,
        milestone
    ):
        self.repository = repository
        self.auth = auth
        self.username = username
        self.state = state
        self.after = after
        self.before = before
        self.milestone = milestone

    def _printReST(self, issues, label):
        """Print section of issues to ``stdout``
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
                raise KeyError(f"Milestone `{self.milestone}` not found")

        return milestone

    def _getDateFromTagOrSHA(self, tagOrSHA):
        """Return date of tag or commit

        If ``tagOrSHA`` is None, return intact to allow all dates.
        If ``tagOrSHA`` corresponds to a ``tag.name``, return date of its commit.
        If ``tagOrSHA`` corresponds to a commit SHA, return date of its commit.
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
        s += [f"{hang}(`#{x.number} <{x.html_url}>`_)"]
        s += ([f"{hang}Thanks to `@{x.user.login} <{x.user.html_url}>`_."]
               if x.user.login not in self.collaborators
               else [])
        return "\n".join(s)

    def format_issue(self, x):
        prefix = "- "
        hang = " " * len(prefix)
        s = [prefix + f"`#{x.number} <{x.html_url}>`_:"]
        s += [textwrap.fill(x.title,
                            initial_indent=hang,
                            subsequent_indent=hang)]
        return "\n".join(s)

    def run(self):
        """Requests issues from GitHub API and prints as reST to ``stdout``
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

def main(
    repository: Annotated[
        str,
        Argument(help="GitHub repository to obtain issues from")
    ] = "usnistgov/fipy",
    token: Annotated[
        str,
        Argument(envvar="FIPY_GITHUB_TOKEN",
                 help="GitHub personal access token with 'repo' scope")
    ] = None,
    username: Annotated[
        str,
        Option(help="GitHub username to authenticate as. "
                    "Superseded by `token`. "
                    "Note: GitHub limits the rate of unauthenticated queries: "
                    "https://developer.github.com/v3/#rate-limiting")
    ] = None,
    state: Annotated[
        str,
        Option(help="Indicates the state of the issues to return. "
                    "Can be either `open`, `closed`, or `all`.")
    ] = "closed",
    after: Annotated[
        str,
        Option(help="Only issues closed at or after this tag, SHA, or date are returned.")
    ] = None,
    before: Annotated[
        str,
        Option(help="Only issues closed at or before this tag, SHA, or date are returned.")
    ] = None,
    milestone: Annotated[
        str,
        Option(help="A string referring to a milestone by its title field. "
                    "If the string `*` is passed, issues with any milestone are accepted. "
                    "If the string `none` is passed, "
                    "issues without milestones are returned. ")
    ] = None,
):
    if token is not None:
        auth = (token,)
    elif username is not None:
        from getpass import getpass

        password = getpass(f"Password for 'https://{self.username}@github.com': ")
        auth = (username, password)
    else:
        raise ValueError("Either a GitHub personal access token or a --username is required")

    return ChangeLog(
        repository=repository,
        auth=auth,
        username=username,
        state=state,
        after=after,
        before=before,
        milestone=milestone
    ).run()

app = typer.Typer()
app.command()(main)

if __name__ == "__main__":
    app()
