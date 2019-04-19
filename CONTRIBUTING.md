## How to contribute to Biotaphy Analyses

#### You found a bug

* Check that it hasn't already be reported by searching our GitHub issues [Issues](https://github.com/biotaphy/analyses/issues).

* If you're unable to find an open issue addressing the problem, [open a new one](https://github.com/biotaphy/analyses/issues/new?assignees=cjgrady&template=bug_report.md). Be sure to include a **title and clear description**, as much relevant information as possible, and a **code sample** or an **executable test case** demonstrating the expected behavior that is not occurring.


#### You wrote a patch for a bug

* Open a new GitHub pull request with the patch.

* Ensure the pull request description clearly describes the problem and solution. Include the relevant issue number if applicable.

* Before submitting, check that your patch follows our coding and testing conventions


#### You want to add a new analysis

* [Submit a new GitHub issue](https://github.com/biotaphy/analyses/issues/new?assignees=&template=feature_request.md) and suggest your analysis.  We want to make sure that it fits before you spend time coding it.

* Write your code and tests following our coding and testing conventions.

* Submit a pull request


#### Coding conventions

* We Use [PEP8](https://www.python.org/dev/peps/pep-0008/)

* Doc strings should follow [Google Style Guidelines](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)

* We use [pytest](https://docs.pytest.org/en/latest/) style tests

* Use pytest with coverage and pep8 to determine if test code is adequate.  At this time, we support Python 2.7 and multiple versions of Python 3.

```
$ py.test tests/ --pep8 analyses -v --cov analyses --cov-report term-missing

$ py.test-3 tests/ --pep8 analyses -v --cov analyses --cov-report term-missing
```

#### You want to update documentation

* Update the documentation in a new branch

* If you update in-line documentation, make sure to rebuild the API doc RST files

```
$ sphinx-apidoc -o ./_sphinx_config/source ./analyses/
```

* If you edit any RST docs (or update API docs), rebuild the html pages

```
$ sphinx-build -b html ./_sphinx_config/ ./docs/sphinx/
```

* Open and submit a new pull request for your updates

Thanks!

BiotaPhy Team
