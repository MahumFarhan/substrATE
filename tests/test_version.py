"""The version number is written in several files; they must agree."""
import os
import re

from click.testing import CliRunner

import substrate
from substrate import cli

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(substrate.__file__)))


def _read(name):
    with open(os.path.join(ROOT, name), encoding='utf-8') as f:
        return f.read()


def _pyproject_version():
    return re.search(r'^version = "([^"]+)"', _read('pyproject.toml'), re.M).group(1)


def test_cli_version_matches_pyproject():
    result = CliRunner().invoke(cli.main, ['--version'])
    assert result.exit_code == 0
    assert result.output.strip().endswith(_pyproject_version())


def test_citation_version_matches_pyproject():
    cff = re.search(r'^version: (\S+)', _read('CITATION.cff'), re.M).group(1)
    assert cff == _pyproject_version()


def test_changelog_has_section_for_this_version():
    assert f'## [{_pyproject_version()}]' in _read('CHANGELOG.md')


def test_patterns_version_matches_pyproject():
    tag = _read(os.path.join('substrate', 'data',
                             'activity_patterns_version.txt')).strip()
    assert tag == f'v{_pyproject_version()}'
