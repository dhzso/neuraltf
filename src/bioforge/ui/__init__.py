"""BioForge UI Layer (Layer 10).

The UI is a thin Streamlit surface sitting on top of the Layer 7 workflow
engine and the Layer 8B evidence framework. The UI itself contains no
omics or evidence logic — it calls :func:`bioforge.workflow.WorkflowExecutor`
and reads the run artifacts the engine writes. This keeps the CLI and UI
behaving identically.

Entry point (after `pip install -e ".[streamlit]"`):

    streamlit run src/bioforge/ui/app.py

The app degrades to a friendly message if the AI provider isn't configured
(StubAssistant is used by Layer 6).
"""
from bioforge.ui.app import main

__all__ = ["main"]
