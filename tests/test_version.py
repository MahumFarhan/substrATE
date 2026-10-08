"""The version number is written in several files; they must agree."""
import os
import re

import pytest
from click.testing import CliRunner

import substrate
from substrate import cli

PKG  = os.path.dirname(os.path.abspath(substrate.__file__))
ROOT = os.path.dirname(PKG)

# The repository files below are only present in a source checkout or an
# editable install, not in a normal installed package.
in_repo = pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, 'pyproject.toml')),
    reason='not running from a source checkout')


def _read(*parts):
    with open(os.path.join(*parts), encoding='utf-8') as f:
        return f.read()


def _cli_version():
    result = CliRunner().invoke(cli.main, ['--version'])
    assert result.exit_code == 0
    return result.output.strip().split()[-1]


def test_cli_reports_a_version():
    assert re.fullmatch(r'\d+\.\d+\.\d+', _cli_version())


def test_patterns_version_matches_cli():
    tag = _read(PKG, 'data', 'activity_patterns_version.txt').strip()
    assert tag == f'v{_cli_version()}'


@in_repo
def test_pyproject_version_matches_cli():
    pyproject = re.search(r'^version = "([^"]+)"',
                          _read(ROOT, 'pyproject.toml'), re.M).group(1)
    assert pyproject == _cli_version()


@in_repo
def test_citation_version_matches_cli():
    cff = re.search(r'^version: (\S+)', _read(ROOT, 'CITATION.cff'), re.M).group(1)
    assert cff == _cli_version()


@in_repo
def test_changelog_has_section_for_this_version():
    assert f'## [{_cli_version()}]' in _read(ROOT, 'CHANGELOG.md')
