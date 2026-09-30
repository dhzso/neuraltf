"""BioForge Streamlit UI.

The UI is a thin Streamlit surface on top of the workflow engine and the
evidence framework. It contains no omics or evidence logic — it calls
:func:`bioforge.workflow.WorkflowExecutor` and reads the run artifacts the
engine writes, which keeps the CLI and UI behaving identically.

Entry point (after `pip install -e ".[streamlit]"`):

    streamlit run src/bioforge/ui/app.py

The app degrades to a friendly message if the AI provider isn't configured
(the assistant falls back to StubAssistant).
"""
from bioforge.ui.app import main

__all__ = ["main"]
