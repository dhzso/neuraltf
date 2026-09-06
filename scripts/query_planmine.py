"""Query PlanMine for functional annotations of NeuralTF candidates.

2026-09-06 audit fix (annotation-universe mismatch): the previous default
input was rank_neural.csv (134 genes) while ALL THREE prioritization
methods score the full rank.csv universe (11,675 genes) — the
GO/domain/ortholog bonuses and the Track-B domain gate were measuring
"was in the old neural subset", not biology (188/230 gate-passing novel
genes had NO annotation rows, and Track-B seats were reshuffled by
annotation POSSESSION). The default query set is now GATE-RELEVANT:

  - all rank_neural.csv genes (neural candidates; keeps prior coverage),
  - PLUS the top-250 novel candidates by base integrated_score (a
    conservative superset of the ~230 gate-passing pool; Track-B seats
    can only ever come from this frontier),
  - PLUS every gene carrying an mmc4 TF flag (the gate's second arm).

Rows for genes already in the parquet are PRESERVED (merged, never
overwritten), and a ``query_status`` base-column separates
``queried_empty`` (no annotations exist) from ``query_failed`` (PlanMine
error — never silently identical to unannotated).

Fetches GO terms, protein domains (Pfam/InterPro), cross-species BLAST
hits and the full transcript sequence for every ``dd_Smed_v6_*``
candidate, and persists:

- ``datasets/processed/planmine_annotations.parquet``  (long-format table)
- ``datasets/processed/planmine_transcripts.fasta``    (sequences, RNAi design)

Usage::

    python scripts/query_planmine.py
    python scripts/query_planmine.py --limit 5     # smoke test
    python scripts/query_planmine.py --out <other.parquet> --fasta <other.fasta)

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


def load_candidates(rank_csv: Path, repo: Path | None = None,
                    n_top_novel: int = 250) -> list[tuple[str, str]]:
    """Build the GATE-RELEVANT query set (2026-09-06 audit fix).

    Members:
      1. every gene in rank_neural.csv (the neural candidates — the
         previous 134-gene universe, fully preserved);
      2. the top ``n_top_novel`` NOVEL candidates (proof_status !=
         known_rnai_validated) by base integrated_score from rank.csv —
         a conservative superset of the Track-B gate-passing pool;
      3. every rank.csv gene carrying a King mmc4 TF flag (the gate's
         second arm) when mmc4 is locatable.

    The set is exactly the genes that can receive a PlanMine-derived
    bonus or gate seat; querying the other ~11k genes is provably inert
    (they can never enter Track B or collect GO/ortholog bonuses).
    """
    df = pd.read_csv(rank_csv)
    if "gene_id" not in df.columns:
        raise ValueError(f"'gene_id' column missing from {rank_csv}")
    name_col = "gene_name" if "gene_name" in df.columns else None

    # 1) neural candidates
    neural_csv = rank_csv.parent / "rank_neural.csv" \
        if rank_csv.name == "rank.csv" else None
    id_set: dict[str, str] = {}
    rows = []
    if neural_csv and neural_csv.exists():
        ndf = pd.read_csv(neural_csv)
        n_name = "gene_name" if "gene_name" in ndf.columns else None
        for _, r in ndf.iterrows():
            gid = str(r["gene_id"]).strip()
            if gid.startswith("dd_Smed_v6_"):
                id_set[gid] = str(r[n_name]).strip() if n_name else ""
    else:
        # called directly on rank_neural.csv: use it as-is
        for _, r in df.iterrows():
            gid = str(r["gene_id"]).strip()
            if gid.startswith("dd_Smed_v6_"):
                id_set[gid] = str(r[name_col]).strip() if name_col else ""
        rows = sorted(id_set.keys())
        return [(g, id_set[g]) for g in rows]

    # 2) top novel frontier by base score
    df["_score"] = pd.to_numeric(df["integrated_score"], errors="coerce").fillna(0)
    novel = df[df["proof_status"] != "known_rnai_validated"] \
        if "proof_status" in df.columns else df
    for _, r in novel.sort_values("_score", ascending=False).head(n_top_novel).iterrows():
        gid = str(r["gene_id"]).strip()
        if gid.startswith("dd_Smed_v6_") and gid not in id_set:
            id_set[gid] = str(r[name_col]).strip() if name_col else ""

    # 3) mmc4 TF-flagged genes (gate's second arm)
    if repo is not None:
        king_dir = repo / "datasets" / "raw" / "Supplementary_Data_ King_2024"
        mmc4 = None
        for cand in (king_dir / "1-s2.0-S2211124724001712-mmc4.xlsx",):
            if cand.exists():
                mmc4 = cand
                break
        if mmc4 is None and king_dir.exists():
            for p in sorted(king_dir.iterdir()):
                if p.suffix.lower() == ".xlsx" and p.stem.lower().endswith("mmc4"):
                    mmc4 = p
                    break
        if mmc4 is not None:
            try:
                tf = pd.read_excel(mmc4, sheet_name="TF")
                flagged = set(tf.loc[tf["TF?"].notna(), "Gene ID"].astype(str))
                name_of = dict(zip(tf["Gene ID"].astype(str),
                                   tf.get("Planarian GenBank Gene Name",
                                          pd.Series(dtype=str)).astype(str)))
                for gid in flagged:
                    gid = gid.strip()
                    if gid.startswith("dd_Smed_v6_") and gid not in id_set:
                        nm = name_of.get(gid, "")
                        id_set[gid] = "" if str(nm) in ("nan", "None") else str(nm)
            except Exception as e:  # noqa: BLE001
                print(f"  (mmc4 TF-flag expansion skipped: {e})")

    return sorted(id_set.items())


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    repo = args.repo.resolve()
    run_dir = repo / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    # 2026-09-06: default = rank.csv (the gate-relevant set is derived
    # from it inside load_candidates); the old default rank_neural.csv
    # covered only 134 of the 11,675-gene scoring universe.
    rank_csv = (args.rank or run_dir / "rank.csv")
    out_parquet = (args.out or repo / "datasets" / "processed" /
                   "planmine_annotations.parquet")
    out_fasta = (args.fasta or repo / "datasets" / "processed" /
                 "planmine_transcripts.fasta")

    candidates = load_candidates(rank_csv, repo=repo)
    if args.limit > 0:
        candidates = candidates[: args.limit]
    print(f"PlanMine annotation query for {len(candidates)} candidates "
          f"(gate-relevant set)")
    print(f"  input : {rank_csv}")
    print(f"  output: {out_parquet}\n          {out_fasta}")

    # Preserve prior coverage: genes already annotated in the existing
    # parquet are NOT re-queried (their rows are kept verbatim); the
    # query set is the gate-relevant set MINUS already-covered genes.
    prior_genes: set[str] = set()
    prior_long: list[dict] = []
    if out_parquet.exists():
        try:
            prior = pd.read_parquet(out_parquet)
            prior_genes = set(prior["gene_id_v6"].astype(str).unique())
            prior_long = prior.to_dict(orient="records")
            print(f"  preserving {len(prior_genes)} already-annotated genes "
                  f"from the existing parquet")
        except Exception as e:  # noqa: BLE001
            print(f"  WARNING: existing parquet unreadable ({e}); "
                  f"starting fresh")
            prior_genes, prior_long = set(), []

    to_query = [(g, n) for g, n in candidates if g not in prior_genes]
    print(f"  to query: {len(to_query)} new genes "
          f"(skipping {len(candidates) - len(to_query)} already covered)")

    client = PlanMineClient(retries=args.retries, rate_limit=args.rate,
                            logger=LOG)

    long_rows: list[dict] = list(prior_long)
    fasta: list[str] = []
    summary: list[dict] = []
    n_missing = 0
    n_failed = 0

    for i, (gid, gname) in enumerate(to_query, start=1):
        query_status = "query_failed"   # set on error below; overwritten
        try:
            ann = client.fetch_contig_annotations(gid)
            query_status = "queried_annotated" if (
                ann["go_terms"] or ann["domains"] or ann["blast_hits"]) \
                else "queried_empty"
        except PlanMineError as exc:
            LOG.error("  [%d/%d] %s FAILED: %s", i, len(to_query), gid, exc)
            # NEVER silently identical to "unannotated": failed queries
            # are flagged so downstream coverage accounting can
            # distinguish them (the 2026-09-04 parquet mixed the two).
            ann = {"contig_id": gid, "length": None, "sequence": None,
                   "go_terms": [], "domains": [], "blast_hits": []}
            n_failed += 1
        seq = ann.get("sequence")
        if not seq:
            n_missing += 1
        # long-format rows: one per annotation, plus a base row
        base = {"gene_id_v6": gid, "gene_name": gname,
                "contig_length": ann["length"]}
        long_rows.append({**base, "kind": "base", "key": "", "value": "",
                          "query_status": query_status})
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
            "query_status": query_status,
        })
        if seq:
            fasta.append(f">{gid} {gname} transcript={gid} length={ann['length']}")
            for j in range(0, len(seq), 60):
                fasta.append(seq[j:j + 60])
        if i % 10 == 0 or i == len(to_query):
            print(f"  [{i}/{len(to_query)}] done "
                  f"(GO/doms/blast, seq_len) last={gid}", flush=True)

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    long_df = pd.DataFrame(long_rows)
    # preserved rows may lack the new column; backfill explicitly
    if "query_status" in long_df.columns:
        long_df["query_status"] = long_df.get("query_status", pd.Series(dtype=str)).fillna("preserved_prior_run")
    long_df.to_parquet(out_parquet, index=False)
    if fasta:
        out_fasta.parent.mkdir(parents=True, exist_ok=True)
        # merge FASTA: keep existing sequences for preserved genes, add new
        existing_txt = out_fasta.read_text(encoding="utf-8") if out_fasta.exists() else ""
        existing_ids = {l.split()[0][1:] for l in existing_txt.splitlines()
                        if l.startswith(">")}
        new_txt = "\n".join(fasta)
        merged = existing_txt.rstrip("\n") + "\n" if existing_txt else ""
        merged += "\n".join(line for line in fasta
                            if not line.startswith(">") or line.split()[0][1:] not in existing_ids)
        out_fasta.write_text(merged + "\n", encoding="utf-8")

    # coverage summary — over the FULL gate-relevant set (queried +
    # preserved), and never counting query_failed as "no annotations"
    dedup = {s["gene_id_v6"]: s for s in summary}
    covered = set(dedup.keys()) | prior_genes
    n_dom_new = sum(1 for s in dedup.values() if s["n_domains"] > 0)
    n_go_new = sum(1 for s in dedup.values() if s["n_go"] > 0)
    print("\n=== PlanMine annotation coverage (gate-relevant set) ===")
    print(f"gate-relevant genes:  {len(candidates)}")
    print(f"  newly queried:     {len(to_query)}")
    print(f"  preserved prior:   {len(prior_genes)}")
    print(f"  covered total:    {len(covered)}")
    print(f"query failures:      {n_failed} (flagged query_status=query_failed)")
    print(f"sequence missing (new): {n_missing}/{len(to_query)}")
    print(f"protein-domain hits (new): {n_dom_new}")
    print(f"GO terms (new):      {n_go_new}")
    print(f"wrote:               {out_parquet}")
    if fasta:
        print(f"wrote:               {out_fasta}")
    return 0


if __name__ == "__main__":
    sys.exit(main())