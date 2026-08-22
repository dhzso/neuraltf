"""Query PlanMine for functional annotations of NeuralTF candidates.

Reads ``rank_neural.csv`` (from a NeuralTF pipeline run), fetches GO terms,
protein domains (Pfam/InterPro), cross-species BLAST hits and the full
transcript sequence for every ``dd_Smed_v6_*`` candidate, and persists:

- ``datasets/processed/planmine_annotations.parquet``  (long-format table)
- ``datasets/processed/planmine_transcripts.fasta``    (sequences, RNAi design)

Usage::

    python scripts/query_planmine.py
    python scripts/query_planmine.py --limit 5     # smoke test
    python scripts/query_planmine.py --out <other.parquet> --fasta <other.fasta>

Logs progress and a coverage summary to stdout.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from bioforge.projects.neuraltf.planmine import PlanMineClient, PlanMineError

LOG = logging.getLogger("query_planmine")
_stderr = logging.StreamHandler(sys.stderr)
_stderr.setFormatter(logging.Formatter("[%(levelname)s] %(name)s: %(message)s"))
LOG.addHandler(_stderr)
LOG.setLevel(logging.INFO)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd(),
                   help="repo root (default: cwd)")
    p.add_argument("--rank", type=Path, default=None,
                   help="path to rank_neural.csv (default: "
                        "projects/NeuralTF/runs/pipeline_run/rank_neural.csv)")
    p.add_argument("--out", type=Path, default=None,
                   help="output parquet (default: datasets/processed/"
                        "planmine_annotations.parquet)")
    p.add_argument("--fasta", type=Path, default=None,
                   help="output FASTA (default: datasets/processed/"
                        "planmine_transcripts.fasta)")
    p.add_argument("--limit", type=int, default=0,
                   help="only query the first N genes (for smoke tests)")
    p.add_argument("--rate", type=float, default=0.25,
                   help="min seconds between requests (rate limiting)")
    p.add_argument("--retries", type=int, default=3)
    return p.parse_args(argv)


def load_candidates(rank_csv: Path) -> list[tuple[str, str]]:
    df = pd.read_csv(rank_csv)
    if "gene_id" not in df.columns:
        raise ValueError(f"'gene_id' column missing from {rank_csv}")
    name_col = "gene_name" if "gene_name" in df.columns else None
    rows = []
    for _, r in df.iterrows():
        gid = str(r["gene_id"]).strip()
        if gid.startswith("dd_Smed_v6_"):
            rows.append((gid, str(r[name_col]).strip() if name_col else ""))
    return rows


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    repo = args.repo.resolve()
    rank_csv = (args.rank or repo / "projects" / "NeuralTF" / "runs" /
                "pipeline_run" / "rank_neural.csv")
    out_parquet = (args.out or repo / "datasets" / "processed" /
                   "planmine_annotations.parquet")
    out_fasta = (args.fasta or repo / "datasets" / "processed" /
                 "planmine_transcripts.fasta")

    candidates = load_candidates(rank_csv)
    if args.limit > 0:
        candidates = candidates[: args.limit]
    print(f"PlanMine annotation query for {len(candidates)} candidates")
    print(f"  input : {rank_csv}")
    print(f"  output: {out_parquet}\n          {out_fasta}")

    client = PlanMineClient(retries=args.retries, rate_limit=args.rate,
                            logger=LOG)

    long_rows: list[dict] = []
    fasta: list[str] = []
    summary: list[dict] = []
    n_missing = 0

    for i, (gid, gname) in enumerate(candidates, start=1):
        try:
            ann = client.fetch_contig_annotations(gid)
        except PlanMineError as exc:
            LOG.error("  [%d/%d] %s FAILED: %s", i, len(candidates), gid, exc)
            ann = {"contig_id": gid, "length": None, "sequence": None,
                   "go_terms": [], "domains": [], "blast_hits": []}
        seq = ann.get("sequence")
        if not seq:
            n_missing += 1
        # long-format rows: one per annotation, plus a base row
        base = {"gene_id_v6": gid, "gene_name": gname,
                "contig_length": ann["length"]}
        long_rows.append({**base, "kind": "base", "key": "", "value": ""})
        for g in ann["go_terms"]:
            long_rows.append({**base, "kind": "go", "key": g["identifier"],
                              "value": g["name"],
                              "namespace": g["namespace"]})
        for d in ann["domains"]:
            long_rows.append({**base, "kind": "domain", "key": d["source"],
                              "value": d["short_name"]})
        for b in ann["blast_hits"]:
            long_rows.append({**base, "kind": "blast", "key": b["target"],
                              "value": b["species"],
                              "description": b["description"]})
        summary.append({
            "gene_id_v6": gid, "gene_name": gname,
            "contig_length": ann["length"],
            "n_go": len(ann["go_terms"]), "n_domains": len(ann["domains"]),
            "n_blast": len(ann["blast_hits"]),
        })
        if seq:
            fasta.append(f">{gid} {gname} transcript={gid} length={ann['length']}")
            for j in range(0, len(seq), 60):
                fasta.append(seq[j:j + 60])
        if i % 10 == 0 or i == len(candidates):
            print(f"  [{i}/{len(candidates)}] done "
                  f"(GO/doms/blast, seq_len) last={gid}", flush=True)

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    long_df = pd.DataFrame(long_rows)
    long_df.to_parquet(out_parquet, index=False)
    if fasta:
        out_fasta.parent.mkdir(parents=True, exist_ok=True)
        out_fasta.write_text("\n".join(fasta) + "\n", encoding="utf-8")

    # coverage summary
    n_hit = len(candidates) - n_missing
    dedup = {s["gene_id_v6"]: s for s in summary}
    n_dom = sum(1 for s in dedup.values() if s["n_domains"] > 0)
    n_go = sum(1 for s in dedup.values() if s["n_go"] > 0)
    n_blast = sum(1 for s in dedup.values() if s["n_blast"] > 0)
    print("\n=== PlanMine annotation coverage ===")
    print(f"candidates:            {len(candidates)}")
    print(f"sequence retrieved:   {n_hit}/{len(candidates)}")
    print(f"protein-domain hits:  {n_dom}/{len(candidates)}")
    print(f"GO terms:             {n_go}/{len(candidates)}")
    print(f"BLAST hits:           {n_blast}/{len(candidates)}")
    print(f"wrote:                {out_parquet}")
    if fasta:
        print(f"wrote:                {out_fasta} ({len(fasta) - len(candidates)} seq lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())