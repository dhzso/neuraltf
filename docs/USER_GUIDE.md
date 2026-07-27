# BioForge User Guide

BioForge is an AI-native bioinformatics workstation built for Deepanshu's
MS thesis (NeuralTF — surfacing high-priority novel neural-fate candidate
transcription factors in planarian single-cell atlases) and reusable for
other cross-atlas TF-prioritization projects.

This guide shows how to use BioForge to its full potential. CLI and UI have
**identical behaviour** — any workflow that runs via the CLI also runs via
the UI; pick whichever surface you prefer.

---

## 1. Quick start (5 minutes)

### 1.1 Start the dev environment

```bash
docker compose up -d
docker exec -w /workspace bioforge-dev bioforge info
```

You should see:

```
bioforge 0.1.0-alpha.1
project: BioForge
...
```

### 1.2 Run the demo pipeline

The demo runs the **entire** evidence→cards→AI pipeline on synthetic data
(no raw datasets required). Both invocations work:

```bash
# CLI
docker exec -w /workspace bioforge-dev bioforge run \
  projects/NeuralTF/workflows/demo_pipeline.yaml \
  --out projects/NeuralTF/runs/demo

# Or, alternatively, edit the workflow and re-run via the UI
# (see §4 below)
```

Output:

```
  ✓ build_demo_records (evidence.demo_rank) in 0.00s
  ✓ build_cards (evidence.build_cards) in 0.00s
  ✓ write_rank (evidence.write_rank_csv) in 0.02s
  ✓ write_cards (report.write_cards_md) in 0.01s
  ✓ summarise (ai.summarize_candidates) in 0.13s
workflow complete: 5 steps run; artifacts in projects/NeuralTF/runs/demo
```

### 1.3 Inspect the artifacts

```bash
docker exec -w /workspace bioforge-dev bash -c "
ls projects/NeuralTF/runs/demo/
# artifacts/  provenance.json  summary.json  ai_summary.md
cat projects/NeuralTF/runs/demo/ai_summary.md
"
```

Each run produces:

| File | Contents |
|------|----------|
| `provenance.json` | Per-step duration, input hashes, STDOUT tail (last 512 chars), loop-item hash if applicable |
| `summary.json`    | High-level run summary (workflow path, step ids, output keys, inputs used) |
| `ai_summary.md`   | Short LLM-authored "what to chase next" paragraph (or stub when AI isn't configured) |
| `artifacts/`      | Empty by default — workflow steps can write outputs here (see §6) |
| `rank.csv` and `evidence_cards.md` | Written by `evidence.write_rank_csv` / `report.write_cards_md` steps into the directory the CLI is run from (or `--out` dir if you update the YAML) |

### 1.4 Launch the UI

```bash
docker exec -w /workspace bioforge-dev streamlit run src/bioforge/ui/app.py
# Open http://localhost:8501
```

The UI has three pages: **Run**, **Results**, **Assistant**.

---

## 2. Architecture at a glance

| Layer | Package | What it does |
|------|---------|--------------|
| 0-5 | `bioforge.core`, `bioforge.plugins`, `bioforge.projects`, `bioforge.cli` | Repo foundation, project scaffolding, CLI shell |
| 6 | `bioforge.ai` | OpenAI-compatible LLM assistant (provider-agnostic) + Stub fallback |
| 7 | `bioforge.workflow` | Declarative YAML workflow engine + step registry + provenance |
| 8A | `bioforge.omics` | Generic scRNA-seq ops (QC, normalize, cluster, trajectory, batch) |
| 8B | `bioforge.evidence` | Cross-atlas evidence framework + cards + confidence tiers |
| 8C | `bioforge.ingest` | Auto-detect dataset ingestion (h5ad / csv / tsv / 10x_mtx / DGE.gz) |
| 9 | `projects/NeuralTF/` | The thesis workflow YAMLs (demo + real-data) |
| 10 | `bioforge.ui` | Streamlit app (Run / Results / Assistant) |

CLI and UI both delegate to Layer 7's `WorkflowExecutor`, so they behave
identically. The UI never imports omics or evidence logic directly.

---

## 3. Running workflows via the CLI

### 3.1 The `bioforge run` command

```bash
bioforge run <workflow.yaml> [--out <dir>] [--input key=value --input key=value ...]
```

- `--out` defaults to `runs/<timestamp>/`. The CLI creates
  `runs/<ts>/provenance.json`, `summary.json`, `ai_summary.md`, and an
  empty `artifacts/` directory.
- `--input key=value` lets you override `$inputs.<key>` references inside
  the workflow on a per-run basis. Repeat the flag for multiple inputs.

### 3.2 Workflow YAML format

```yaml
description: My workflow

inputs: {disease_classifier: neural}

steps:
  - id: pre
    target: qc
    inputs: {adata: $inputs.x}
    params: {max_pct_mt: 15.0}

  - id: cluster
    target: cluster
    inputs: {adata: $pre.adata}     # $step.output reference
    params: {resolution: 0.6}

  - id: per_atlas
    target: evidence.score_per_atlas
    loop: $inputs.atlas_list        # iterate over this list (literal or $ref)
    inputs: {adata: $item}          # $item resolves to current loop item
```

Reference forms:
- `{"$step": "pre", "$output": "adata"}` (dict form)
- `"$pre.adata"` (string form)
- `"$inputs.<key>"` (workflow-level input override)
- `"$item"` (current loop iteration value)

### 3.3 Registered steps (Layer 7 + 8A + 8B + 8C + AI + workflow.steps)

To list every available step at runtime:

```python
from bioforge.workflow.registry import StepRegistry
print(StepRegistry.instance().known())
```

The names you'll see:

- **8A omics**: `ingest`, `qc`, `normalize`, `cluster`, `trajectory`,
  `batch_correct`, `storage.save_anndata`
- **8C ingest**: `ingest.fincher`, `ingest.plass`, `ingest.king`
- **8B evidence**: `evidence.demo_rank`, `evidence.write_rank_csv`,
  `evidence.load_bridge_csv`, `evidence.score_per_atlas`,
  `evidence.combine_scores`, `evidence.add_rnai_stream`,
  `evidence.add_correlation_stream`, `evidence.add_function_stream`,
  `evidence.assign_tiers`, `evidence.build_cards`
- **Reports**: `report.write_cards_md`
- **AI**: `ai.summarize_candidates`

### 3.4 CLI subcommands (parity with UI features)

```bash
bioforge info                      # Version + resolved config
bioforge datasets list [--category raw|processed|reference|cache]
bioforge projects list             # List project scaffolds
bioforge projects create <name>   # Scaffold a new project (charter layout)
bioforge plugins list              # Entry-points registered as plugins
bioforge run <wf.yaml> [...]       # Execute a YAML workflow
```

---

## 4. Running workflows via the Streamlit UI

1. Launch the UI (see §1.4).
2. **Run page**: pick a workflow YAML from `projects/<project>/workflows/*`.
   Paste one dataset source per line in the text area. A **format badge**
   appears next to each line telling you whether the file is a recognised
   `.h5ad`, DGE `.txt.gz`, 10x mtx directory, etc. Lines of the form
   `fincher=path/to/file` become `$inputs.fincher` references in the
   workflow; plain lines become `$inputs.dataset_<N>`.
3. Click **Run workflow**. A live log of step completions appears. On
   success the latest run directory is remembered and the **Results** page
   becomes populated.
4. **Results page**: pick a run from the dropdown; the rank CSV renders as
   a sortable table. A dropdown lets you **inspect one candidate** which
   expands the corresponding evidence-card markdown fragment inline.
5. **Assistant page**: chat with the configured LLM. If you type
   `tool: <name> arg1|arg2` the UI dispatches the named tool locally and
   prints its result. Available tools: `lookup_gene`, `summarize_candidates`,
   `inspect_anndata`. When no API key is set, the `StubAssistant` returns
   deterministic placeholder replies so the page is never broken.

---

## 5. Configuring the AI provider

BioForge's AI Layer is **provider-agnostic via the OpenAI-compatible HTTP
API** — works with OpenAI, OpenRouter, Together, Groq, Anyscale, vLLM,
SGLang, LM Studio, or Ollama's OpenAI shim. Swap provider by changing the
base URL.

```bash
docker compose down
export BIOFORGE_AI_API_KEY=sk-...
export BIOFORGE_AI_BASE_URL=https://api.openai.com/v1  # default
export BIOFORGE_AI_MODEL=gpt-4o-mini                  # default
docker compose up -d
```

To use, say, **OpenRouter**:

```bash
export BIOFORGE_AI_BASE_URL=https://openrouter.ai/api/v1
export BIOFORGE_AI_MODEL=anthropic/claude-3.5-sonnet
export BIOFORGE_AI_API_KEY=or-...
```

To use **Ollama** (local, free):

```bash
export BIOFORGE_AI_BASE_URL=http://host.docker.internal:11434/v1
export BIOFORGE_AI_MODEL=llama3
export BIOFORGE_AI_API_KEY=dummy-anything     # Ollama doesn't check it
```

If no API key is exported, BioForge transparently falls back to
`StubAssistant` and all workflows still complete end-to-end — this is
your thesis reproducibility guarantee: reviewers can run your workflow
without an API key.

---

## 6. Running the real-data NeuralTF pipeline

Once you have the Fincher / Plass / King supplementary files locally:

```bash
docker exec -w /workspace bioforge-dev bioforge run \
  projects/NeuralTF/workflows/real.yaml \
  --out projects/NeuralTF/runs/real_$(date +%Y%m%d-%H%M%S) \
  --input fincher=datasets/raw/GSE111764_GEO_Fincher_atlas/GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz \
  --input plass=datasets/processed/plass.h5ad \
  --input king_catalog=datasets/raw/Supplementary_Data_\ King_2024/1-s2.0-S2211124724001712-mmc4.xlsx \
  --input king_rnai=datasets/raw/Supplementary_Data_\ King_2024/1-s2.0-S2211124724001712-mmc5.xlsx \
  --input king_corr=datasets/raw/Supplementary_Data_\ King_2024/1-s2.0-S2211124724001712-mmc6.xlsx \
  --input king_atlas=datasets/raw/Supplementary_Data_\ King_2024/1-s2.0-S2211124724001712-mmc7.xlsx \
  --input bridge=projects/NeuralTF/data/bridge.csv
```

The real pipeline steps (in order):

1. **ingest_fincher**, **ingest_plass**, **ingest_king** — Load DGE
   matrices + King supplementary xlsx tables (`tf_catalog`, `rnai`,
   `correlations`, `atlas`).
2. **qc** — Layer 8A QC pipeline (filter + QC metrics). The same step can
   be looped across Atlases via the `loop:` keyword.
3. **load_bridge** — Load a `v4↔v6 gene bridge` CSV. (Build the CSV
   separately via `build_bridge_from_names` — see §7.).
4. **score_per_atlas** — Loop over atlases: per-cluster log2FC +
   Shapiro-based specificity score per TF (restricted to King's 716
   candidate TF catalogue when available).
5. **combine_scores** — Aggregate per-atlas records by gene ID, applying
   `REPRODUCIBILITY = n_atlases_supported / 3` per candidate.
6. **add_rnai_stream** — Score `RNai = 1.0` if the gene appears in the King
   mmc5 RNAi table, `0.0` otherwise.
7. **add_correlation_stream** — Score `CORRELATION = (g0_corr − x1_corr) ·
   3, capped to 1.0` per King mmc6. (The "maturation" gain Deepanshu's
   thesis cares about.)
8. **add_function_stream** — Score `FUNCTION` from
   `ontology.annotate_function` (soxB / pou4l-1 / myoD → categorised;
   unknown → 0).
9. **assign_tiers** — Run `confidence.assign_tiers` for the integrated
   score + supporting-streams count → high / medium / low.
10. **build_cards** — Build per-candidate `EvidenceCard`s with:
    - `integrated_score`, `tier`, `proof_status`
    - `supporting_streams`
    - `per_source` breakdown with notes (e.g. which atlases a TF was
      enriched in)
    - `atlases_supported` (e.g. `["fincher", "plass", "king"]`)
    - `suggested_followups` (wet-lab / ortholog suggestions)
11. **write_rank** + **write_cards** — Persist CSV + markdown.
12. **summarise** — `ai.summarize_candidates` produces the one-paragraph AI
    summary saved to `ai_summary.md`.

### 6.1 Proof status classes (thesis-critical)

Every evidence card carries a `ProofStatus`:

- **`known_rnai_validated`** — existing RNAi phenotype in King mmc5. Action:
  reuse for combinatorial-code studies.
- **`novel_candidate`** — no RNAi, no prior FSTF record. **This is the
  class Deepanshu's wet-lab experimentation follows up on.**
- **`prior_fstf_not_tested`** — known FSTF in literature, but King's assay
  did not show a phenotype. Worth re-testing with different markers / cell
  types.

The three classes share the same confidence score; the **separation
appears in the markdown card and the UI**, not in the integrated score —
so a really strong novel candidate can sit at the top of the list ahead
of weaker RNAi-validated known TFs.

---

## 7. Building the v4↔v6 gene bridge CSV

BioForge never guesses gene ID equivalents from numeric prefixes. You
build an explicit bridge table:

```python
from bioforge.evidence import build_bridge_from_names
import pandas as pd

v6 = pd.DataFrame({
    "gene_name": ["soxB", "myoD", "pitx", "pou4l-1"],
    "v6_id":     ["dd_Smed_v6_15104_0_1", "dd_Smed_v6_22001_0_1",
                  "dd_Smed_v6_x",        "dd_Smed_v6_y"],
})
v4 = pd.DataFrame({
    "gene_name": ["soxB", "myoD", "pitx", "hnf4"],
    "v4_id":     ["dd_Smed_v4_6001", "dd_Smed_v4_2", "dd_Smed_v4_3", "dd_Smed_v4_9001"],
})

bridge = build_bridge_from_names(v6, v4)
bridge.df.to_csv("projects/NeuralTF/data/bridge.csv", index=False)
print(f"bridged {bridge.n_bridged}/{bridge.n_rows} names")
```

Output formats:

- `bridge.csv` with columns `gene_name, v6_id, v4_id` — some `v4_id` may
  be `NaN` for genes only present in the v6 build.

For real-scale bridging, the `gene_name` column on both sides should come
from the supplementary gene-name tables shipped by Fincher (`GSE111764_*/family.soft.gz`)
and King (`mmc2.xlsx`); both are loaded by the `bioforge.evidence.readers`
package.

---

## 8. Umplug-and-go: running without an AI key

```bash
unset BIOFORGE_AI_API_KEY
docker compose up -d
docker exec -w /workspace bioforge-dev bioforge run projects/NeuralTF/workflows/demo_pipeline.yaml --out /tmp/somerun
```

The workflow still runs end-to-end; the `ai.summarize_candidates` step
yields the `StubAssistant`'s deterministic placeholder (`[bioforge-stub]
AI provider not configured. Would have answered: "..."`). Your code never
hard-crashes on AI being offline — this is intentional for thesis
reproducibility.

---

## 9. Dataset ingestion

### 9.1 Auto-detected formats

```python
from bioforge.ingest import ingest_dataset

adata = ingest_dataset("/home/deepanshu/data/x.h5ad")
adata = ingest_dataset("datasets/raw/.../dge.txt.gz")
adata = ingest_dataset("https://example.org/atlas.csv")
adata = ingest_dataset("GSE12345")  # see §9.2 below
```

BioForge sniffs file content (HDF5 magic, gzip magic, text-tab vs text-comma)
rather than relying on extensions only.

### 9.2 Accessions

Direct GEO/SRA fetch isn't wired into the framework yet — the workflow's
`ingest` step raises a friendly `UnknownFormatError` advising you to fetch
the data first via `geo` (NCBI Gene Expression Omnibus `wget` helper) or
`sra-toolkit fasterq-dump`. BioForge then loads the local file via
auto-detection. This keeps BioForge out of the rapidly-changing
NCBI-protocol business.

### 9.3 FASTQ → matrix (optional)

A stub `bioforge.ingest.fastq` declares the future recipe; raising
`NotImplementedError` until the optional `[fastq]` extra (kb-python OR
salmon) plus recipe wiring lands. For now, run FASTQ → 10x_mtx externally
and pass the resulting directory directly to `ingest_dataset`.

---

## 10. Programmatic API examples

### 10.1 Score one TF manually

```python
from bioforge.evidence import (
    EvidenceRecord, EvidenceSource, EvidenceScorer,
    assign_tiers, build_evidence_card, render_card_markdown,
)

r = EvidenceRecord(gene_id="dd_Smed_v6_15104_0_1", gene_name="soxB")
r.add_score(EvidenceSource.EXPRESSION,      1.0,  note="log2FC=5.07,atlas=king")
r.add_score(EvidenceSource.SPECIFICITY,     0.8,  note="entropy=0.22")
r.add_score(EvidenceSource.REPRODUCIBILITY, 1.0,  note="atlases=3/3")
r.add_score(EvidenceSource.RNai,            0.0,  note="not in mmc5")
r.add_score(EvidenceSource.CORRELATION,     0.6,  note="g0_corr=0.41")
r.add_score(EvidenceSource.FUNCTION,        1.0,  note="category=neural")

card = build_evidence_card(r)
print(card.tier.value)         # 'high'
print(card.proof_status.value) # 'novel_candidate'
print(render_card_markdown(card))
```

### 10.2 Custom workflow step

```python
from bioforge.workflow.registry import register

@register("my.print_genes",
          inputs=["records"],
          outputs=["printed"],
          description="Print top-10 candidate TFs.")
def print_genes(records, top_n: int = 10) -> dict:
    for r in records[:top_n]:
        print(f"{r.gene_name}: {r.supporting_streams()} streams")
    return {"printed": top_n}
```

Now use `target: my.print_genes` in any YAML workflow.

### 10.3 AI tool dispatch programmatically

```python
from bioforge.ai import inspect_anndata, lookup_gene, summarize_candidates
import json

print(json.loads(inspect_anndata("runs/latest/x_processed.h5ad")))
print(json.loads(lookup_gene("soxB", bridge_path="projects/NeuralTF/data/bridge.csv")))
print(json.loads(summarize_candidates("runs/latest/rank.csv", top_n=5)))
```

---

## 11. Reproducibility guarantees

Every `bioforge run` invocation writes:

- `provenance.json` — step id, target, started_at, duration, input_hash,
  stdout_tail (last 512 chars), loop_item hash, error.
- `summary.json` — overall run description (workflow YAML, step ids,
  per-run `--input` overrides).
- `ai_summary.md` — LLM-authored summary (or stub).

Anyone can reproduce your run by:

1. Checking out the same commit.
2. Using the same `docker/base/requirements.lock.txt` (172 pinned deps).
3. Running the same `bioforge run <yaml> --input ...` command found in
   `summary.json`.

The two key reproducibility rules the framework enforces for you:
- The bridge table is **explicit** (you supply `--input bridge=...`); no
  numeric-prefix guessing.
- AI is **opt-in**: the StubAssistant is the default if no key is set, so
  the workflow always completes.

---

## 12. Testing your own hypotheses / extending BioForge

### 12.1 Add a new analysis step

```bash
# 1. Write the callable
echo '
from bioforge.workflow.registry import register

@register("custom.filter_to_known_neural",
          inputs=["records"], outputs=["records"],
          description="Keep only TFs whose function category is 'neural'.")
def fn(records: list) -> dict:
    from bioforge.evidence.schema import EvidenceSource
    kept = [r for r in records if r.notes.get(EvidenceSource.FUNCTION, "").endswith("=neural")]
    return {"records": kept}
' > src/bioforge/workflow/custom_steps.py
```

```yaml
# 2. Use it in workflow YAML
steps:
  - id: filter_neural
    target: custom.filter_to_known_neural
    inputs: {records: $assign_tiers.records}
    params: {}
```

### 12.2 Add a new AI tool

```python
from bioforge.ai.tools import register_tool

@register_tool("count_neuronal_candidates")
def count_neuronal_candidates(rank_csv: str) -> str:
    import json, csv
    n = sum(1 for row in csv.DictReader(open(rank_csv)) if "neural" in row.get("gene_name",""))
    return json.dumps({"n": n})
```

Use it from the Assistant page in the UI as `tool: count_neuronal_candidates
runs/latest/rank.csv`.

### 12.3 Run the test suite yourself

```bash
docker exec -w /workspace bioforge-dev python -m pytest tests/ -q
# 182 passed, 2 warnings in 40.34s
```

The warnings are upstream pandas deprecation notices in scanpy
(`Pandas4Warning: The copy keyword is deprecated…`) and don't affect any
test outcome.

---

## 13. Troubleshooting

### "AI provider not configured"

This is the StubAssistant working as designed — see §5 to configure a
real provider.

### `UnknownFormatError` from `ingest_dataset`

The function surfaced gracefully rather than crashing. Two usual causes:
- The file is a `matrix.mtx` outside its directory — collect it together
  with `barcodes.tsv` + `features.tsv` in one directory.
- The accession (`GSE12345`) — see §9.2.

### Workflow hangs forever

Each step in the executor runs synchronously. With `progress_cb`,
the CLI/UI shows you which step is currently running. If a step runs
longer than expected (e.g. CellRank on a large atlas), check the last
`INFO` log line — the executor has no built-in timeouts; you'd add a
Python-level timeout in the step itself if desired.

### `BridgeTable missing required columns`

The bridge CSV must have `gene_name, v6_id, v4_id`. Re-read §7.

### UI doesn't show my new run

The Results page looks under `projects/NeuralTF/runs/`. If you ran with
`--out /tmp/xyz`, that directory is outside the UI's view. Either pass an
`--out` under `projects/NeuralTF/runs/` or load the markdown file
manually.

---

## 14. Where to look next

- `docs/decisions/0001-bioforge-layer-plan-v2.md` — overall layered plan
- `docs/decisions/0002-layer-8b-evidence-integration-design.md` — Evidence
  framework module design (weights, tiers, bridge-table contract)
- `docs/decisions/0003-ai-workflow-ui-design.md` — AI Layer + workflow
  engine + UI design
- `CHANGELOG.md` — every notable change since the repo's start
- `projects/NeuralTF/README.md` — end-to-end thesis project manifest

Good luck with the wet-lab validation. The novel candidate TFs surfaced
under `proof_status = novel_candidate` at tier `high` should be your first
RNAi + FISH targets.
