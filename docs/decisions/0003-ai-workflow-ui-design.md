# ADR-0003: AI Layer, Workflow Engine, Dataset Ingestion, UI, and Cross-Validation Design

- **Status:** Accepted
- **Date:** 2026-07-27
- **Deciders:** Deepanshu (Owner), OpenCode (acting Architecture Lead)
- **Relates to:** ADR-0001 (layer plan), ADR-0002 (8B evidence framework)

## Context

After Layers 0–5 and 8A/8B we have:
- A working Docker dev environment and CLI shell.
- Reusable omics ops (`bioforge.omics.{qc,normalize,cluster,trajectory,batch}`).
- An Evidence Integration Framework (`bioforge.evidence.*`).

What's still missing for the thesis is **the experience end-to-end**: a
researcher types an SRA accession, the workstation ingests the data, runs
the omics pipeline, integrates cross-atlas evidence, and returns a ranked
list of neural-fate candidate TFs — clearly separating *novel untested* from
*RNAi-validated known* candidates, feeding both classes back into the same
confidence scheme. The user (per Q&A) wants the AI assistant built into the
loop, an OpenAI-compatible provider model, a Streamlit-style interactive UI
(not a desktop GUI), and robust dataset ingestion that won't choke on
unexpected file types.

## Decisions

### Layer 6 — AI Layer (`bioforge.ai`)

- **Provider-agnostic via OpenAI-compatible HTTP**: the consumer of any
  plugin/provider only ever calls `AIAssistant.complete(messages=...)` and
  gets back a `ChatResponse`. The first concrete backend is an
  `OpenAICompatClient` that POSTs to `${base_url}/v1/chat/completions` with
  an `Authorization: Bearer ${api_key}` header. This speaks the protocol
  used by OpenAI, OpenRouter, Together, Groq, Anyscale, Mistral,
  vLLM/SGLang, Ollama (with the OpenAI shim), and LM Studio — so a
  researcher swaps providers by editing `~/.config/bioforge/ai.toml`.
- **Configuration discovery**: read from env vars first
  (`BIOFORGE_AI_BASE_URL`, `BIOFORGE_AI_API_KEY`, `BIOFORGE_AI_MODEL`),
  then a YAML/TOML config file, then a built-in default base_url that
  points to OpenAI's public endpoint.
- **Stub fallback**: if no key is configured, a `StubAssistant` returns
  canned responses so tests and offline development don't fail. Calling
  `complete()` raises `AIProviderNotConfiguredError` only when an explicit
  "must use live" flag is set on the request — by default we degrade to
  the stub so the workflow never hard-fails on AI being unavailable.
- **Assistant-facing tools**: a small set of deterministic tools are wired
  up so the LLM can call them via a structured function-call protocol:
  - `inspect_anndata(path)` — shape, obs/var columns, layer names.
  - `summarize_ranked_candidates(rank_csv)` — returns a one-paragraph
    summary of the top-K candidate TFs.
  - `lookup_gene(name_or_id)` — looks up a gene in the BridgeTable /
    TF catalog.
  These tools are Python callables registered with the AI module so the
  assistant is genuinely useful to a researcher, not a blind chatbot.

### Layer 7 — Workflow Engine (`bioforge.workflow`)

- **Declarative YAML files** describe a workflow as an ordered list of
  named steps; each step names a target (a registered callable) plus
  `inputs` (refs to upstream step outputs or workflow-level inputs) and
  `params` (kwargs). Outputs are bound to names so downstream steps can
  reference them (`step_id.output_name`).
- **`StepRegistry`**: callables register with `@register("qc")`-style
  decorator. Each registered callable declares its **input contract** (an
  AnnData path, a DataFrame, etc.) via Python type hints — the executor
  uses these hints to coerce inputs (e.g. unwrap `DataRef`s to actual
  values).
- **Executor**: runs steps in dependency-resolved order; captures
  `ExecutionContext` per step (duration, STDOUT lines, error trace,
  hashable input fingerprint); writes a `provenance.json` per run that
  records every step's target version, params, and input hash.
- **CLI**: `bioforge run workflow.yaml [--input SRA123 --out runs/]` invokes
  the executor, streams progress to STDERR, and writes results to
  `runs/<timestamp>/{artifacts, provenance.json, ai_summary.md}`.
- **Composition**: a workflow can reference another workflow as a sub-step
  so Layer 9's NeuralTF pipeline is just a workflow consuming 8A/8B.

### Layer 8C — Dataset Ingestion (`bioforge.ingest`)

- **Single entry point**: `ingest_dataset(source, *, dest_dir=None) -> AnnData`.
  `source` may be a geographic/SRA accession (GSE/GSM/SRP/SRR), a URL, or
  a local path.
- **Format auto-detection** by content sniffing **not extension only**:
  | Format                             | Signature                              | Reader                       |
  |------------------------------------|----------------------------------------|------------------------------|
  | gzipped DGE txt (Fincher-style)    | gzip → first bytes are tab-separated   | `read_dge_gz`                |
  | 10x matrix market (h5 + mtx)       | directory with `matrix.mtx` + `barcodes.tsv` | `scanpy.read_10x_mtx` |
  | `.h5ad`                            | HDF5 magic (`\x89HDF\r\n\x1a\n`)       | `anndata.read_h5ad`          |
  | CSV matrix (genes × cells)         | first line is comma-sep header         | `pandas.read_csv` → AnnData  |
  | TSV matrix                         | first line is tab-sep header           | `pandas.read_csv(sep="\t")`  |
- **Unreadable inputs degrade gracefully**: `UnknownFormatError` is raised
  only when *no* reader accepts the file; the workflow catches it and emits
  a warning + a placeholder AnnData with `uns['ingestion_error']` set so
  the run doesn't crash. The user-visible path (UI) catches this exception
  and renders a clear "couldn't identify this dataset format, please upload
  one of: …" message instead of a stack trace.
- **Optional FASTQ→matrix orchestration**: a sub-module `bioforge.ingest.fastq`
  declares a Snakemake-style *recipe* that the workflow can opt into via
  `recipe: fastq_to_10x`. The recipe knows about `kb-python` (kallisto |
  bustools) for 10x-style and `salmon` for non-UMI bulk. Actual execution
  of that recipe requires the optional `[fastq]` extra and writes a 10x-mtx
  directory, which the auto-detector then reads normally. The first
  iteration is a thin orchestrator around the recommended toolchain — not
  performance-tuned — but it gives a clear runnable path.

### Layer 8B extension — Cross-Atlas Cross-Validation + Novel vs Validated

Extend `bioforge.evidence.confidence` with:

- A dedicated `ProofStatus` enum: `KNOWN_RNAI_VALIDATED`, `NOVEL_CANDIDATE`,
  `PRIOR_FSTF_NOT_TESTED`. New helper classifies each candidate based on
  the RNAi stream score and the TF-catalog prior FSTF flag.
- `EvidenceCard` dataclass that carries: `RecordRef`, `Tier`, `IntegratedScore`,
  `ProofStatus`, `atlases_supported: set[str]` (e.g. `{"fincher","plass","king"}`),
  `per_source_breakdown: dict[EvidenceSource, tuple[float, str]]` (score + note),
  `suggested_followups: list[str]` (suggested next experiments).
- `build_evidence_card(record, scorer, bridge, atlases, rnai_table)` produces
  one card. `build_cards_for_records(...)` produces a list.
- A `render_card_markdown(card)` helper that writes a single markdown
  fragment per candidate — used both for `runs/.../ai_summary.md` and
  for the Streamlit UI evidence panel.
- The thesis-critical rule: **novel untested candidates are surfaced
  separately** so a researcher can pick window FYIs, while the RNAi-validated
  stream contributes to the integrated score exactly as before (no junk
  signal discarded).

### Layer 9 — NeuralTF workflow (`projects/NeuralTF/`)

- A new `projects/NeuralTF/` directory is scaffolded by `bioforge projects
  create NeuralTF` (existing Layer 5 capability).
- The project commits `workflows/discover_neural_tfs.yaml` referencing
  the same step ids shipped in Layers 7/8. Steps:
  1. `ingest.fincher`, `ingest.plass`, `ingest.king`
  2. `qc`, `normalize`, `cluster`, `trajectory` per AnnData (looped by the
     workflow engine's `loop:` keyword — see Layer 7)
  3. `evidence.load_bridge_csv` → `evidence.score_per_atlas` →
     `evidence.combine_scores`
  4. `evidence.assign_tiers`
  5. `evidence.build_cards` → `report.write_cards_md` → optionally
     `ai.summarize_candidates`
- ais integrated into the workflow via the same `step` interface; the LLM
  gets the ranked candidate list and produces a single-paragraph "what to
  chase next" note saved alongside the CSV.

### Layer 10 — Streamlit UI (`bioforge.ui`)

- One file: `src/bioforge/ui/app.py` — registered as a `[streamlit]` extra.
- Three pages via sidebar: **Run**, **Results**, **Assistant**.
  - **Run**: a text box accepting SRA/GSE/GSM accessions or local file
    paths; "Load datasets" button; for each dataset shows a detected
    format badge; "Run workflow" button invokes the Layer 7 executor
    with a `ProgressCallback` that updates a Streamlit progress bar.
  - **Results**: reads the latest run dir; renders a sortable candidate
    table (gene, name, tier, integrated score, proof status); clicking
    a row expands to show the `EvidenceCard` markdown and an optional
    AI summary section.
  - **Assistant**: chat panel wired to `AIAssistant.complete`; greets
    with a stub message when not configured; supports the registered
    tools so the user can ask "show me the top 5 neural TFs and their
    evidence" and the model returns tool-call requests the UI executes.
- The UI is purely an orchestration surface — it never imports omics or
  evidence logic directly. All operations go through `bioforge.workflow` /
  `bioforge.evidence` so the CLI and UI behave identically.

## Alternatives Considered

1. **LiteLLM instead of raw OpenAI HTTP** — gives a roster of providers out
   of the box but adds a non-trivial dependency we don't otherwise need.
   The OpenAI compat approach gives us 90% of providers with one HTTP
   client and zero new deps.

2. **Nextflow instead of custom workflow YAML** — too heavy for a thesis
   project; we'd be locked into a runtime that doesn't fit our embedded
   Python container well. Our tiny YAML + Python-callables registry is
   enough for what NeuralTF needs and trivially extensible.

3. **No AI Layer 6 (defer)** — explicit user request flagged AI integration
   as a thesis selling point; building the base layer now keeps the
   open-source ambition real without forcing every workflow to use LLMs.

## Consequences

- Adds runtime deps: `requests` (Layer 6), `pyyaml` is already present,
  `streamlit` (Layer 10) under a new `[streamlit]` extra so the headless
  install footprint stays minimal.
- The UI keeps a declarative boundary: any callable registered with Layer
  7 is automatically usable from both CLI and Streamlit; we add new
  workflow steps as plain Python functions and they appear in the UI
  without UI changes.
- We make AI **optional and graceful**: workflows that don't need an LLM
  still run end-to-end with zero API configuration, which matters for
  reproducible thesis artifacts.
