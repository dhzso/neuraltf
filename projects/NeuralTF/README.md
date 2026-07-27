# NeuralTF — Planarian Neural-Fate Transcription Factor Discovery

The First Consumer for BioForge. Goal: surface high-confidence novel
(neuron-fate specifying) transcription factor candidates from the
Fincher 2018 + Plass 2018 + King 2024 planarian cell atlases, ranked by
multi-source evidence so Deepanshu's wet-lab follow-up (RNAi + FISH)
prioritises the most defensible hits.

## Layout

```
projects/NeuralTF/
  README.md                  # this file
  workflows/
    discover_neural_tfs.yaml # the rank-and-evidence pipeline
  runs/                       # workflow run artifacts (gitignored)
```

## Running the pipeline

The workflow is declarative YAML understood by `bioforge run` (Layer 7):

```bash
bioforge run projects/NeuralTF/workflows/discover_neural_tfs.yaml \
  --out projects/NeuralTF/runs/$(date +%Y%m%d-%H%M%S)
```

The demo run uses the `evidence.demo_rank` synthetic candidate generator
so the end-to-end pipeline works without raw data:

```bash
bioforge run projects/NeuralTF/workflows/demo_pipeline.yaml \
  --out projects/NeuralTF/runs/demo
```

## Real-data pipeline (under construction)

The full pipeline ingests three raw datasets, runs QC/normalize/cluster/
PAGA on each, loads the King TF catalog + RNAi + correlation tables,
builds an explicit v4↔v6 bridge, scores every candidate TF across six
evidence streams, ranks by integrated score, and emits per-candidate
evidence cards separating novel candidates from RNAi-validated known FSTFs.

See `workflows/real.yaml` once the bridge CSV is pinned down; building
the bridge from Fincher/King name-matching is a project step tracked
separately.

## AI assistance

Workflow steps can call `ai.summarize_candidates`; this uses the
Provider-AI Layer 6 to ask an LLM for a one-paragraph "what to chase
next" note. With no API key configured the StubAssistant returns a
deterministic placeholder so the workflow never hard-fails.
