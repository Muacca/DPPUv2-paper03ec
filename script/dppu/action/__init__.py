"""
DPPU Action Layer
=================

Lagrangian, effective potential, stability analysis, and EC action builder.
"""

from .lagrangian import compute_lagrangian, compute_action
from .potential import compute_effective_potential, get_potential_function, subs_zero_modes
from .stability import analyze_stability, find_equilibrium_r, scan_vacuum_3d, StabilityType
from .ec_action import compute_c2_ec, build_veff_ec, build_veff_ec_func

__all__ = [
    'compute_lagrangian',
    'compute_action',
    'compute_effective_potential',
    'get_potential_function',
    'subs_zero_modes',
    'analyze_stability',
    'find_equilibrium_r',
    'scan_vacuum_3d',
    'StabilityType',
    'compute_c2_ec',
    'build_veff_ec',
    'build_veff_ec_func',
]
