"""Heuristics for classifying a user-supplied dataset source string."""
from __future__ import annotations

import re

ACCESSION_PATTERNS = (
    r"^GSE\d+$",
    r"^GSM\d+$",
    r"^SRP\d+$",
    r"^SRR\d+$",
    r"^SRX\d+$",
    r"^ERP\d+$",
    r"^DRP\d+$",
)


def is_accession(source: str) -> bool:
    return any(re.match(p, source.strip()) for p in ACCESSION_PATTERNS)


def is_url(source: str) -> bool:
    return source.startswith(("http://", "https://", "ftp://"))


def is_local_path(source: str) -> bool:
    # Heuristic: contains a slash, backslash, or has a path-like extension.
    return any(c in source for c in ("/", "\\")) or source.startswith(".")


def classify_source(source: str) -> str:
    """Return one of 'accession', 'url', 'local_path'."""
    if is_url(source):
        return "url"
    if is_accession(source):
        return "accession"
    return "local_path"
