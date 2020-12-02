Releasing a SMUTHI version
===========================

1. Run tests
------------

Use nosetests to automatically run all tests on your local computer. This step is not strictly necessary, because automatic testing is also covered by the GitLab CI pipeline after every commit.

2. Update version number
------------------------

The place to set the version number is smuthi/version.py as well as docs/conf.py

3. Update requirements
----------------------

If a new third party package was added, on which SMUTHI now depends: Add that requirement in setup.py

4. Update README.rst
--------------------

Add a new section "What's new in version xyz"

Update the list of people who have contributed.

5. Update docs/about_smuthi.rst
-------------------------------

Update the list of people who have contributed.

6. Push to the online repository
--------------------------------

(Or create a merge request)

7. Add a release tag in the online repository
---------------------------------------------

This will trigger the following actions:

- The GitLab CI pipeline will create a source distribution as well as Linux binary wheels and upload everything to PyPi
- The appveyor CI pipeline will create Windows binary wheels and upload them to PyPi
