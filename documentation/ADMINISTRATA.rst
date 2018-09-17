=============
Git practices
=============

Refer to :ref:`documentation:GIT` for the current branching conventions.

--------
Branches
--------

Whether fixing a bug or adding a feature, all work on FiPy should be based
on a reported `GitHub issue`_. Assuming issue number 12345, branch the code::

    $ BRANCH=issue12345-Summary_of_what_branch_addresses
    $ git checkout -b $BRANCH develop

Edit and add to branch::

    $ emacs ...
    $ git commit -m "refactoring_stage_A"
    $ emacs ...
    $ git commit -m "refactoring_stage_B"

Merging changes from develop to the branch
------------------------------------------

Make sure ``develop`` is up to date::

    $ git fetch origin

Merge updated state of ``develop`` to the branch::

    $ git diff develop
    $ git merge develop

Resolve any conflicts and test::

    $ python setup.py test

Submit branch for code review
-----------------------------

.. attention::

   **Administrators Only!**

   Push the code to GitHub for automated testing::

       $ git push origin $BRANCH

   Check the Buildbot_ status. Fix (or, if absolutely necessary, document)
   any failures.

Paste the result of::

    $ git request-pull develop origin

(or whatever is appropriate for your branch and clone) into a comment on
the ticket this branch addresses. You can send a message to the mailing
list about it if you like, but the FiPy developers should see the pull
request via RSS feed.

.. note::

   A ``request-pull`` requires that your repository be publicly accessible.
   If that's not possible, then alternatively, you can do a::

       $ git format-patch

   and then attach the :file:`*.patch` files to the ticket.

Refactoring complete: merge branch to develop
---------------------------------------------

.. attention::

   **Administrators Only!**

   First, follow the instructions for
   `Merging changes from develop to the branch`_.

   Merge the branch to ``develop``::

       $ git checkout develop
       $ git diff $BRANCH
       $ git merge $BRANCH

   Resolve any conflicts and test::

       $ python setup.py test

   Push the code to GitHub for automated testing::

       $ git push origin develop

   When completely done with the branch::

       $ git branch -D $BRANCH
       $ git push origin :$BRANCH

---------
Bug fixes
---------

.. note::

   This is not well thought out. By and large, we don't do this.

At the point some fix is made to an old version n.m *that is not at the tip
of ``master``*, make a branch from that old release (this step not
necessary if the branch already exists due to a previous fix)::

    $ git branch version-n_m refs/tags/version-n_m

Proceed as with other Branches_, but instead branching from ``develop``,
do development work off of the historical branch::

    $ BRANCH=ticket12345-Summary_of_what_branch_addresses
    $ git checkout -b $BRANCH version-n_m

Edit and commit as usual.

If appropriate, after successful code review and merger to the
``version-n_m`` branch, the changes should also be merged to ``develop``::

    $ git checkout develop
    $ git merge version-n_m

.. attention::

   When complete, the ``version-n_m`` branch is not merged to ``master``.

================
Making a Release
================

.. attention::

   **Administrators Only!**

------
Source
------

Make sure ``develop`` is ready for release::

   $ git checkout develop

Check items in the issues_ and update the :file:`README.txt`::

   $ git commit README.txt -m "REL: update new features for release"

.. attention::

   If Buildbot_ doesn't show all green boxes for this release, make sure to
   add appropriate notes in :file:`README.txt` or :file:`INSTALLATION.txt`!

-------------------
Release from master
-------------------

::

    $ git checkout master
    $ git merge develop

Resolve any conflicts and push to ``master``::

    $ git tag --annotate version-x_y master
    $ git push --tags origin master

Clean the working copy::

    $ git clean -fd

.. note::

   Alternatively, clone into a clean repository.

Build the documentation and the web pages::

    $ python setup.py bdist_egg
    $ python setup.py build_docs --pdf --html --cathartic

Build the compressed distribution::

    $ rm MANIFEST
    $ python setup.py sdist

Test the installed compressed distribution::

    $ cpvirtualenv trunk test
    $ mkdir tmp
    $ cd tmp
    $ cp ../dist/FiPy-${FIPY_VERSION}.tar.gz .
    $ tar zxvf FiPy-${FIPY_VERSION}.tar.gz
    $ cd FiPy-${FIPY_VERSION}
    $ python setup.py install
    $ cd ..
    $ python -c "import fipy; fipy.test()"
    $ deactivate
    $ rmvirtualenv test
    $ cd ..
    $ \rm -rf tmp

-------
Windows
-------

Build a windows executable installer::

    $ rm MANIFEST
    $ python setup.py bdist --formats=wininst

Combine the windows installer and examples into one archive::

    $ rm MANIFEST
    $ FIPY_VERSION=XXX
    $ ln dist/FiPy-${FIPY_VERSION}.win32.exe .
    $ cp MANIFEST.in MANIFEST.in.bkup
    $ cp MANIFEST-WINDOWS.in MANIFEST.in
    $ python setup.py sdist --dist-dir=dist-windows --formats=zip
    $ cp MANIFEST.in.bkup MANIFEST.in
    $ unlink FiPy-${FIPY_VERSION}.win32.exe
    $ mv dist-windows/FiPy-${FIPY_VERSION}.zip dist/FiPy-${FIPY_VERSION}.win32.zip

------
Debian
------

Make sure stdeb_ and debhelper_ are installed::

    $ cd CLEAN
    $ python setup.py --command-packages=stdeb.command bdist_deb
    $ mv deb_dist/python-fipy_${FIPY_VERSION}-1_all.deb dist/python-fipy_${FIPY_VERSION}-1_all.deb

------
Upload
------

Tag the repository as appropriate (see `Git practices`_ above).

Upload the build products to PyPI

    $ python setup.py sdist upload

Upload the build products and documentation from :file:`dist/` and
the web site to CTCMS ::

    $ export FIPY_WWWHOST=bunter:/u/WWW/wd15/fipy
    $ export FIPY_WWWACTIVATE=updatewww
    $ python setup.py upload_products --pdf --html --tarball --winzip

.. warning:: Some versions of ``rsync`` on Mac OS X have caused problems
   when they try to upload erroneous ``\rsrc`` directories. Version 2.6.2
   does not have this problem.

Make an announcement to `fipy@nist.gov`_

Build (``python setup.py bdist --formats=wininst``) a Windows `PyVTK`_
executable and upload to download page.

.. _GitHub issue: https://github.com/usnistgov/fipy/issues/new
.. _issues: https://github.com/usnistgov/fipy/issues
.. _Buildbot: http://build.cmi.kent.edu:8010/tgrid
.. _fipy@nist.gov: mailto:fipy@nist.gov
.. _PyVTK: http://cens.ioc.ee/projects/pyvtk/
.. _stdeb: http://github.com/astraw/stdeb
.. _debhelper: http://kitenet.net/~joey/code/debhelper/
