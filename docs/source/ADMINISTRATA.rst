=============
Git practices
=============

Refer to :ref:`documentation:GIT` for the current branching conventions.

--------
Branches
--------

Whether fixing a bug or adding a feature, all work on FiPy should be
conducted on a branch and submitted as a `pull request`_. If there is
already a reported GitHub_ issue_, name the branch accordingly::

    $ BRANCH=issue12345-Summary_of_what_branch_addresses
    $ git checkout -b $BRANCH master

Edit and add to branch::

    $ emacs ...
    $ git commit -m "refactoring_stage_A"
    $ emacs ...
    $ git commit -m "refactoring_stage_B"

Merging changes from master to the branch
-----------------------------------------

Make sure ``master`` is up to date::

    $ git fetch origin

Merge updated state of ``master`` to the branch::

    $ git diff origin/master
    $ git merge origin/master

Resolve any conflicts and test::

    $ python setup.py test

Submit branch for code review
-----------------------------

If necessary, fork_ the `fipy repository`_.

Add a "remote" link to your fork::

    $ git remote add <MYFORK> <MYFORKURL>

Push the code to your fork on GitHub_::

    $ git push <MYFORK> $BRANCH

Now `create a pull request`_ from your ``$BRANCH`` against the ``master``
branch of ``usnistgov/fipy``.  The `pull request`_ should initiate
automated testing.  Check the `Continuous Integration`_ status.  Fix (or,
if absolutely necessary, document) any failures.

.. note::

   If your branch is still in an experimental state, but you would like to
   check its impact on the tests, you may prepend "``WIP:``" to your `pull
   request`_ title.  This will prevent your branch from being merged before
   it's complete, but will allow the automated tests to run.

   Please be respectful of the `Continuous Integration`_ resources and do
   the bulk of your testing on your local machine or against your own
   `Continuous Integration`_ accounts (if you have a lot of testing to do,
   before you create a `pull request`_, push your branch to your own
   fork_ and enable the `Continuous Integration`_ services there.

You can avoid testing individual commits by adding "``[skip ci]``" to the
commit message title.

When your `pull request`_ is ready and successfully passes the tests, you
can `request a pull request review`_ or send a message to the mailing list
about it if you like, but the FiPy developers should automatically see the
pull request and respond to it without further action on your part.

Refactoring complete: merge branch to master
--------------------------------------------

.. attention::

   **Administrators Only!**

Use the GitHub_ interface to `merge the pull request`_.

.. note::

   Particularly for branches with a long development history, consider
   doing a `Squash and merge`_.


.. _CONTINUOUSINTEGRATION:

======================
Continuous Integration
======================

| |Tests|_ |Documentation|_ |nix|_

We use the :term:`Azure` and :term:`GitHub Actions` cloud services for
:term:`Continuous Integration` (CI).  These CIs are configured in
|.azure/pipelines.yml|_, |.github/workflows/NISTtheDocs2Death.yml|_, and
|.github/workflows/nix.yml|_.

.. note::

   In order to focus on breakages introduced by changes to :term:`FiPy`, a
   `pull request`_ is normally built with one of the :ref:`CONDALOCKFILES`,
   whereas the nightly builds use an :file:`environment.yml` in order to
   catch breakages introduced by :term:`FiPy`'s prerequisites.

   A `pull request`_ may be tested with the latest prerequisites by setting
   the ``CONDA_ENVIRONMENT_NOT_LOCK`` environment variable in
   `Azure at queue time`_.

.. |Tests|         image:: https://dev.azure.com/guyer/FiPy/_apis/build/status/usnistgov.fipy?branchName=master
.. _Tests:         https://dev.azure.com/guyer/FiPy/_build?definitionId=2
.. |Documentation| image:: https://github.com/usnistgov/fipy/actions/workflows/NISTtheDocs2Death.yml/badge.svg
.. _Documentation: https://github.com/usnistgov/fipy/actions/workflows/NISTtheDocs2Death.yml
.. |nix|           image:: https://github.com/usnistgov/fipy/actions/workflows/nix.yml/badge.svg
.. _nix:           https://github.com/usnistgov/fipy/actions/workflows/nix.yml

.. |.azure/pipelines.yml| replace::    :file:`{FiPySource}/.azure/pipelines.yml`
.. _.azure/pipelines.yml: https://github.com/usnistgov/fipy/blob/master/.azure/pipelines.yml
.. |.github/workflows/NISTtheDocs2Death.yml| replace::    :file:`{FiPySource}/.github/workflows/NISTtheDocs2Death.yml`
.. _.github/workflows/NISTtheDocs2Death.yml: https://github.com/usnistgov/fipy/blob/master/.github/workflows/NISTtheDocs2Death.yml
.. |.github/workflows/nix.yml| replace::    :file:`{FiPySource}/.github/workflows/nix.yml`
.. _.github/workflows/nix.yml: https://github.com/usnistgov/fipy/blob/master/.github/workflows/nix.yml
.. _Azure at queue time: https://learn.microsoft.com/en-us/azure/devops/pipelines/process/variables?view=azure-devops&tabs=yaml%2Cbatch#allow-at-queue-time

.. _CONDALOCKFILES:

===============
Conda Lockfiles
===============

The `conda-lock <https://github.com/conda/conda-lock>`_ lockfiles in
:file:`environments/locks/` can be updated with::

    $ for solver in petsc pysparse scipy trilinos
      do
        conda-lock lock \
          --file environments/${solver}-environment.yml \
          --lockfile environments/locks/conda-${solver}-lock.yml
        conda-lock render \
          --filename-template environments/locks/conda-${solver}-{platform}.lock \
          environments/locks/conda-${solver}-lock.yml
      done

.. attention::

   Do not merge new lockfiles to ``master`` without validating that
   everything still works.

================
Making a Release
================

.. attention::

   **Administrators Only!**

------
Source
------

Make sure ``master`` is ready for release::

   $ git checkout master

Check the issue_ list and update the :ref:`CHANGELOG`::

   $ git commit CHANGELOG.txt -m "REL: update new features for release"

.. note::

   You can use::

      $ python setup.py changelog --after=<x.y>

   or::

      $ python setup.py changelog --milestone=<x.z>

   to obtain a ReST-formatted list of every GitHub_ `pull request`_ and issue_
   closed since the last release.

   Particularly for major and feature releases, be sure to curate the
   output so that it's clear what's a big deal about this release.
   Sometimes a `pull request`_ will be redundant to an issue_, e.g.,
   "``Issue123 blah blah``".  If the `pull request`_ fixes a bug,
   preference is given to the corresponding issue_ under **Fixes**.
   Alternatively, if the `pull request`_ adds a new feature, preference is
   given to the item under **Pulls** and corresponding issue_ should be
   removed from **Fixes**.  If appropriate, be sure to move the "Thanks to
   @mention" to the appropriate issue_ to recognize outside contributors.

   ..  attention:: Requires PyGithub_ and Pandas_.

.. attention::

   If `Continuous Integration`_ doesn't show all green boxes for this
   release, make sure to add appropriate notes in :file:`README.txt` or
   :file:`INSTALLATION.txt`!

.. _PyGithub: https://pygithub.readthedocs.io
.. _Pandas: https://pandas.pydata.org

-------------------
Release from master
-------------------

::

    $ git checkout master

Resolve any conflicts and tag the release as appropriate (see `Git
practices`_ above)::

    $ git tag --annotate x.y master

Push the tag to GitHub_::

    $ git push --tags origin master

Upon successful completion of the `Continuous Integration`_ systems, fetch
the tagged build products from Azure_ Artifacts and place in
:file:`{FiPySource}/dist/`:

 * :file:`dist-Linux/FiPy-{x.y}-none-any.whl`
 * :file:`dist-Linux/FiPy-{x.y}.tar.gz`
 * :file:`dist-Windows_NT/FiPy-{x.y}.zip`
 * :file:`dist-docs/FiPy-{x.y}.pdf`
 * :file:`dist-docs/html-{x.y}.tar.gz`

From the :file:`{FiPySource}` directory, unpack :file:`dist/html-{x.y}.tar.gz`
into :file:`docs/build` with::

    $ tar -xzf dist/html-{x.y}.tar.gz -C docs/build

.. _Azure:         https://dev.azure.com/guyer/FiPy/_build?definitionId=2

------
Upload
------

Attach
 * :file:`dist/FiPy-{x.y}-none-any.whl`
 * :file:`dist/FiPy-{x.y}.tar.gz`
 * :file:`dist/FiPy-{x.y}.zip`
 * :file:`dist/FiPy-{x.y}.pdf`

to a `GitHub release`_ associated with tag x.y.

Upload the build products to PyPI with twine_::

    $ twine upload dist/FiPy-${FIPY_VERSION}.*

Upload the web site to CTCMS ::

    $ export FIPY_WWWHOST=bunter:/u/WWW/wd15/fipy
    $ export FIPY_WWWACTIVATE=updatewww
    $ python setup.py upload_products --html

.. warning:: Some versions of ``rsync`` on Mac OS X have caused problems
   when they try to upload erroneous ``\rsrc`` directories. Version 2.6.2
   does not have this problem.

.. _GitHub release: https://github.com/usnistgov/fipy/releases

----------------------------
Update conda-forge feedstock
----------------------------

Once you push the tag to GitHub_, the fipy-feedstock_ should automatically
receive a pull request.  Review and amend this pull request as necessary
and ask the `feedstock maintainers`_ to merge it.

This automated process only runs once an hour, so if you don't wish to wait
(or it doesn't trigger for some reason), you can manually generate a pull
request to update the fipy-feedstock_ with:

* revised version number
* revised sha256 (use ``openssl dgst -sha256 /path/to/fipy-x.y.tar.gz``)
* reset build number to ``0``

--------
Announce
--------

Make an announcement to `fipy@list.nist.gov`_

.. _GitHub: https://github.com/
.. _fipy repository: https://github.com/usnistgov/fipy
.. _issue: https://github.com/usnistgov/fipy/issues
.. _pull request: https://github.com/usnistgov/fipy/pulls
.. _fork: https://help.github.com/en/articles/fork-a-repo
.. _create a pull request: https://help.github.com/en/articles/creating-a-pull-request
.. _request a pull request review: https://help.github.com/en/articles/requesting-a-pull-request-review
.. _merge the pull request: https://help.github.com/en/articles/merging-a-pull-request
.. _Squash and merge: https://help.github.com/en/articles/about-pull-request-merges/#squash-and-merge-your-pull-request-commits
.. _twine: https://pypi.org/project/twine
.. _fipy-feedstock: https://github.com/conda-forge/fipy-feedstock
.. _fipy@list.nist.gov: mailto:fipy@list.nist.gov
.. _feedstock maintainers: https://github.com/conda-forge/fipy-feedstock#feedstock-maintainers
