# Getting Started

**ssw-py** is designed to be simple to install and use.

```{contents}
```

## Motivation

`ssw-py` exists to help oligo alignment in python.  The upstream [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) is quality cross-platform solution to use as the foundation of this tool.

## Requirements

**ssw-py** is built and tested on MacOS, Linux and Windows 64-bit systems; we do not provide official Windows support. Python versions 3.8 - 3.11 builds are supported.

Wheels are released for CPython versions following the [EOL model](https://devguide.python.org/versions/).


## Installation

If you want to install the latest stable build of **ssw-py**, you can
install it using `pip`:

```bash
$ pip install ssw-py
```

**NOTE**: We support wheel builds for PyPi for the 3 most recent CPython versions. Target platforms for wheels are MacOS `x86-64` `arm64`, Linux `x86-64`, and Windows `x86-64`.

If your Python version and platform fall outside this such as Linux `aarch64` it is confirmed `ssw-py` builds on this platform but it is not supported as our build GitHub actions runners do not run these builds expediently.

## Alignment

The {py:mod}`ssw` module includes support for

### Workflow

The easiest way to run the code is
Example usage:

```python
```

## Advanced Installation

Users interested in contributing to development may want to work with the
latest development build. To get the latest and greatest code, head over
[our Github repo](https://github.com/libnano/ssw-py) and clone the
repo or download a tarball. Building from source is easy.

If you don't install the latest build via pip or conda, you might have to install
`Cython`, prior to running the `setup.py` script:

```
$ pip install Cython
```

Or via `conda`:

```
$ conda install Cython
```

Then run:

```
$ python setup.py install
```

or if you are developing `ssw-py` enhancements:

```
$ python setup.py build_ext --inplace
```

We recommend running `setup.py` with either `build_ext --inplace` or
`develop` rather than `install` if you are testing development builds.
`build_ext --inplace` will build the Cython and C API extensions in the
package directory without copying any files to your local environment
site-packages directory (so you can import and run tests from within the
package) and `develop` will build in place and then put symlinks in your
site packages directory.


## Testing

Every commit pushed to
[the ssw-py GitHub repo](https://github.com/libnano/ssw-py) is tested to
ensure it builds properly and passes our unit testing framework as a GitHub action

If you'd like to run the tests yourself, we suggest the following workflow:

```
$ git clone https://github.com/libnano/ssw-py
$ cd ssw-py
$ python setup.py build_ext --inplace
$ pytest
```

NOTE: `pip` / `conda` install `pytest` if not in your environment


## Contributing

Contributions are welcomed via pull requests.

Contributions that improve the stability, compatibility or test coverage are
best way to interface with the project and dev team.  Feature requests via
GitHub issues are also welcome.

Contact the `ssw-py` maintainers prior to beginning your work to make sure
it makes sense for the project.

By contributing, you also agree to release your code under the MIT License

After a successful PR will be listed under the [contributors](https://github.com/libnano/ssw-py/graphs/contributors).


### Forking

A forking workflow is preferred for all pull requests.

### Branch naming

Branch naming is preferred to use the format:

```
<GitHub user-name>-<short keyword description of change>
```

Keep branch names not too long.  A good example would be for the user `grinner`
for a documentation update for the 1.0.0 staging branch:

```bash
$ git checkout -b grinner-docs-update-1.0.0-pass-01
```

With the trailing 01 indicative of it being part of several potential

Another example pass that focuses on code clarity comments would be:

```bash
$ git checkout -b grinner-code-clarity-and-comments
```

### Development

Development requires the use of C Python 3.8+, [pytest](https://docs.pytest.org) and
[pre-commit](https://pre-commit.com) as they are used to build and run `ssw-py`
code CI in the GitHub Action.

Install these dependencies in your python development environment
(`virtualenv`, `conda`, etc):

```bash
$ pip install cython pre-commit pytest
# or
$ conda install cython pre-commit pytest
```

Install `pre-commit` in repo the with:

```bash
$ pre-commit install
```

To ensure the git hook is excecuted on every commit.

### Pull Requests

Pull Requests should meet the following requirements:

1. Excellent PR description describing all changes made. Please use markdown syntax highlighting to help readability.
2. If change is code related, have test coverage for the changes implemented.
3. Attempt to make the PR 1 commit only. Multiple are OK if it helps illustrate the change better.
4. Commit messages should describe the changes.
5. Provided you contact the maintainers in advance, theu will code review your PR, provide feedback and squash merge your code on approval.

**TIP**: Interactive `rebase` is helpful to fix old commit messages.
For example, run:

```bash
$ git rebase -i HEAD~2
```

To rebase the last 2 commits. Use `s` to mark the most recent commit(s), save, then
modify the collective commit messages to update poor commit messages.
