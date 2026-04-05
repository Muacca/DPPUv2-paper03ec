"""
DPPU Utils: Shared Utility Functions
====================================

Provides common utility functions used across the package.

Modules:
- levi_civita: Levi-Civita symbols (3D and 4D)
- symbolic: Symbolic computation helpers (prove_zero, witness search)
- block_classifier: Curvature block classification
- parameters: Unified parameter interface
- printing: Console output helpers (hline, print_header, print_sub)
"""

from .levi_civita import epsilon_symbol, levi_civita_4d
from .symbolic import prove_zero, find_nonzero_witness, generate_test_points
from .printing import hline, print_header, print_sub

__all__ = [
    'epsilon_symbol',
    'levi_civita_4d',
    'prove_zero',
    'find_nonzero_witness',
    'generate_test_points',
    'hline',
    'print_header',
    'print_sub',
]
