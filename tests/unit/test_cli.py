"""Smoke tests for the bioforge CLI via Click's CliRunner."""
from click.testing import CliRunner

from bioforge import __version__
from bioforge.cli.main import cli


def test_version_option() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_help_lists_subcommands() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "info" in result.output
    assert "datasets" in result.output
    assert "projects" in result.output
    assert "plugins" in result.output


def test_info_command(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["info"])
    assert result.exit_code == 0
    assert "bioforge" in result.output
    assert "project:" in result.output


def test_datasets_list_on_empty_repo(tmp_path, monkeypatch) -> None:
    # Set up a minimal dataset tree
    (tmp_path / "datasets" / "raw").mkdir(parents=True)
    (tmp_path / "datasets" / "processed").mkdir(parents=True)
    (tmp_path / "datasets" / "reference").mkdir(parents=True)
    (tmp_path / "datasets" / "cache").mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["datasets", "list"])
    assert result.exit_code == 0
    assert "raw/" in result.output


def test_projects_list_empty(tmp_path, monkeypatch) -> None:
    (tmp_path / "projects").mkdir()
    monkeypatch.chdir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["projects", "list"])
    assert result.exit_code == 0
    assert "no projects" in result.output


def test_projects_create(tmp_path, monkeypatch) -> None:
    (tmp_path / "projects").mkdir()
    monkeypatch.chdir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["projects", "create", "TestProj"])
    assert result.exit_code == 0
    assert "Created project 'TestProj'" in result.output
    assert (tmp_path / "projects" / "TestProj" / "data").is_dir()


def test_plugins_list_no_plugins(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["plugins", "list"])
    assert result.exit_code == 0
    assert "no plugins registered" in result.output or "bioforge" not in result.output.split("\n", 1)[0]
