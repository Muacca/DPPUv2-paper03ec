"""Notebook helper: load and display the saved manuscript figures."""
from pathlib import Path
from paper_results.figures import NAMES,figure

DEFAULT_DATA=Path(__file__).resolve().parents[1]/'data/paper03ec_v4'

def show(name,source=DEFAULT_DATA):
    """Return a matplotlib figure. This function never recomputes the model."""
    return figure(name,source)
