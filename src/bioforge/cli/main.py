"""BioForge CLI entry point.

Built with Click. The CLI is intentionally minimal at Layer 5 — it exposes
project plumbing (config, datasets, projects, plugins) but no analysis
subcommands (those arrive with Layer 8 and beyond).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click

from bioforge import __version__
from bioforge.core.config import load_config
from bioforge.core.datasets import DatasetManager
from bioforge.core.logging import configure_logging
from bioforge.projects.manager import ProjectManager


def _repo_root() -> Path:
    """Return the repo root inferred from the current working directory."""
    cwd = Path.cwd()
    return cwd  # BioForge is invoked from the repo root


@click.group()
@click.version_option(__version__, prog_name="bioforge")
@click.option(
    "--config",
    "config_path",
    type=click.Path(dir_okay=False),
    default=None,
    help="Path to a YAML configuration file.",
)
@click.pass_context
def cli(ctx: click.Context, config_path: str | None) -> None:
    """BioForge — AI-native bioinformatics workstation."""
    cfg = load_config(config_path)
    configure_logging(cfg.logging)
    ctx.ensure_object(dict)
    ctx.obj["config"] = cfg
    ctx.obj["repo_root"] = _repo_root()


@cli.command()
@click.pass_context
def info(ctx: click.Context) -> None:
    """Show BioForge version and resolved configuration."""
    cfg = ctx.obj["config"]
    click.echo(f"bioforge {__version__}")
    click.echo(f"project: {cfg.project}")
    click.echo(f"logging.level: {cfg.logging.level}")
    click.echo(f"logging.file: {cfg.logging.file}")
    click.echo(f"datasets.root: {cfg.datasets.root}")


@cli.group()
def datasets() -> None:
    """Inspect the dataset layout."""


@datasets.command("list")
@click.option(
    "--category",
    type=click.Choice(["raw", "processed", "reference", "cache"]),
    default=None,
    help="List only one category.",
)
@click.pass_context
def datasets_list(ctx: click.Context, category: str | None) -> None:
    """List datasets in each (or one) category."""
    cfg = ctx.obj["config"]
    root = ctx.obj["repo_root"]
    mgr = DatasetManager(root, cfg.datasets)
    cats = [category] if category else ["raw", "processed", "reference", "cache"]
    for cat in cats:
        items = mgr.list(cat)
        click.echo(f"{cat}/ ({len(items)} items)")
        for name in items:
            click.echo(f"  - {name}")


@cli.group()
def projects() -> None:
    """Inspect or create research projects."""


@projects.command("list")
@click.pass_context
def projects_list(ctx: click.Context) -> None:
    """List available research projects."""
    mgr = ProjectManager(ctx.obj["repo_root"])
    names = mgr.list()
    if not names:
        click.echo("(no projects)")
        return
    for name in names:
        try:
            p = mgr.resolve(name)
            click.echo(f"- {p.name}  ({p.root})")
        except Exception as exc:  # noqa: BLE001 - report, don't crash listing
            click.echo(f"- {name}  (invalid: {exc})")


@projects.command("create")
@click.argument("name")
@click.pass_context
def projects_create(ctx: click.Context, name: str) -> None:
    """Create a new project scaffold."""
    mgr = ProjectManager(ctx.obj["repo_root"])
    try:
        p = mgr.create(name)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)
    click.echo(f"Created project '{p.name}' at {p.root}")


@cli.group()
def plugins() -> None:
    """Inspect registered plugins."""


@plugins.command("list")
@click.pass_context
def plugins_list(ctx: click.Context) -> None:
    """List discovered plugins (entry points only, does not load them)."""
    from bioforge.plugins import PluginManager

    mgr = PluginManager()
    eps = mgr.discover()
    if not eps:
        click.echo("(no plugins registered)")
        return
    for ep in eps:
        click.echo(f"- {ep.name}  ({ep.value})")


@cli.command("run")
@click.argument(
    "workflow_yaml",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--out",
    "out_dir",
    type=click.Path(),
    default=None,
    help="Output directory for run artifacts (default runs/<timestamp>).",
)
@click.option(
    "--input",
    "inputs",
    multiple=True,
    help=(
        "Per-run input override as key=value. Repeat for multiple inputs "
        "(e.g. --input fincher=path/dge.txt.gz --input plass=path.h5ad). "
        "Accessible inside the workflow as $inputs.fincher."
    ),
)
@click.pass_context
def run_cmd(ctx: click.Context, workflow_yaml: str, out_dir: str | None,
            inputs: tuple[str, ...]) -> None:
    """Execute a declarative YAML workflow end-to-end."""
    import datetime as dt
    import json
    import logging

    # Lazy import so omics deps don't have to be importable for `bioforge info`.
    from bioforge.workflow import WorkflowExecutor, WorkflowRun
    import bioforge.workflow.steps as _steps  # noqa: F401 — registers steps

    ts = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    out_path = Path(out_dir or f"runs/{ts}")
    artifacts_dir = out_path / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    log = logging.getLogger("bioforge.cli.run")

    # Parse --input key=value overrides
    extra_inputs: dict[str, Any] = {}
    for kv in inputs:
        if "=" not in kv:
            raise click.UsageError(f"--input must be 'key=value'; got: {kv}")
        k, v = kv.split("=", 1)
        extra_inputs[k] = v
    if extra_inputs:
        log.info("per-run inputs: %s", extra_inputs)

    run = WorkflowRun.from_yaml(workflow_yaml)
    log.info("loaded workflow with %d steps from %s", len(run.steps), workflow_yaml)

    def progress(step_id: str, target: str, duration: float) -> None:
        click.echo(f"  ✓ {step_id} ({target}) in {duration:.2f}s")

    executor = WorkflowExecutor(progress_cb=progress)
    try:
        outputs = executor.execute(run, extra_inputs=extra_inputs)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"workflow failed: {exc}", err=True)
        sys.exit(1)

    provenance_p = out_path / "provenance.json"
    provenance_p.write_text(json.dumps(executor.provenance, indent=2), encoding="utf-8")
    log.info("wrote provenance to %s", provenance_p)

    # AI summary (if the workflow exposed 'summary' on the ai.summarize step)
    summary_text = None
    for step_id, out in outputs.items():
        if isinstance(out, dict) and "summary" in out:
            summary_text = out["summary"]
            break
    if summary_text:
        (out_path / "ai_summary.md").write_text(f"# AI summary\n\n{summary_text}",
                                                encoding="utf-8")

    summary = {
        "workflow_yaml": str(workflow_yaml),
        "out_dir": str(out_path),
        "n_steps": len(run.steps),
        "step_ids": [s.id for s in run.steps],
        "outputs": {k: list(v.keys()) for k, v in outputs.items()},
        "inputs": extra_inputs,
    }
    (out_path / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    click.echo(f"workflow complete: {len(outputs)} steps run; artifacts in {out_path}")


def main() -> None:
    """Entry point used by the ``bioforge`` console script."""
    cli()  # pylint: disable=no-value-for-parameter


if __name__ == "__main__":
    main()
