"""Dataset readers for the evidence-integration framework.

Each reader module loads one source's supplementary files into tidy
:class:`pandas.DataFrame` objects so the rest of the framework is
source-agnostic.
"""
from bioforge.evidence.readers import king, fincher, plass

__all__ = ["king", "fincher", "plass"]
