# =============================================================================
# DPPUv2 Interactive Phase Diagram & Potential Viewer (v4ec — EC-Weyl Engine)
# =============================================================================
# EC extension of the v4 viewer.  UI is identical; only the engine backend
# and the Weyl-status indicator are updated for Einstein-Cartan theory.
#
# Prerequisites:
#   - %matplotlib widget  (run before importing this module)
#   - dppu package on sys.path  (add script/ two levels up)
#
# EC changes vs v4 (LC):
#   - Engine: UnifiedEngine (dppu.topology.unified) replaces topology-specific
#     S3S1Engine / T3S1Engine / Nil3S1Engine from _DPPUv2_Phase2/
#   - Nil³ radial symbol: params['R'] (scale param); S3/T3 use params['r']
#   - T³: alpha not a free symbol (flat → C²=0 always); alpha ignored silently
#   - get_weyl_c2_status: Nil³ is NON-conformally flat even at ε=0
#     (C²_LC = 4/(3R⁴) ≠ 0 — EC paper Theorem 1 / Sec. 5.1)
#
# Features (unchanged from v4):
#   - Topology: S³×S¹ / T³×S¹ / Nil³×S¹
#   - Torsion Mode: MX / VT / AX  (on-demand, cached)
#   - NY Variant: FULL / TT / REE
#   - ε (squashing): geometric deformation of S³ and Nil³
#   - α_W (Weyl coupling): adds α·C² term to Lagrangian
#   - Weyl status indicator: dynamic C²=0 / C²≠0 feedback
#   - Phase diagram (V, η) at fixed (θ_NY, ε, α_W)
#   - Potential V_eff(r) comparison for up to 3 selected points
#   - Click-to-select points on phase diagram
#   - Axis scale switching (linear / log / symlog)
#
# UI Layout:
#   ┌─ Configuration ──────────────────────────────────┐  ┌─ Status ─────┐
#   │ [Topology ▼] [Mode ▼] [NY Variant ▼]             │  │ status msg   │
#   │ θ_NY: ━━━━━━━━━━━━━━━━━━━━━━━━━ [1.0]           │  │ cache info   │
#   │ ┌─ Geometric Parameters ─────────────────────┐  │  └──────────────┘
#   │ │ ε (squashing): ━━━━━━━━━━━━━━━━━ [0.00]    │  │
#   │ │ α_W (Weyl):    ━━━━━━━━━━━━━━━━━ [0.00]    │  │
#   │ │ [Weyl status indicator]                     │  │
#   │ └────────────────────────────────────────────┘  │
#   └──────────────────────────────────────────────────┘
#   ┌─ Phase Diagram Controls ─────┐ ┌─ Potential Controls ──────────┐
#   │ V max / η min / η max        │ │ r max / |V_eff| range         │
#   │ Point selection (3 points)   │ │ Axis scale toggles            │
#   │ [Draw Phase Diagram]         │ │ [Update Potential]            │
#   └──────────────────────────────┘ └──────────────────────────────┘
#   ┌─ Output (matplotlib figure) ──────────────────────────────────────┐
#   │  [Phase Diagram (V, η)]         [Effective Potential V_eff(r)]   │
#   └───────────────────────────────────────────────────────────────────┘
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import ipywidgets as widgets
from ipywidgets import Layout, HBox, VBox, Output
from IPython.display import display, clear_output
from scipy.optimize import minimize_scalar
import sys
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Path Setup: ensure dppu package (in _DPPUv2_Phase2/) is importable
# =============================================================================

def _setup_dppu_path():
    """Add _DPPUv2_Phase2/ to sys.path so that `import dppu` works."""
    try:
        # Works when this file is imported as a module (__file__ is set)
        _here = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        # Fallback: assume the notebook's cwd is scripts/visualize/
        _here = os.getcwd()
    _root = os.path.abspath(os.path.join(_here, '..', '..'))
    if _root not in sys.path:
        sys.path.insert(0, _root)

_setup_dppu_path()

# =============================================================================
# Engine Imports (v4 dppu package)
# =============================================================================

try:
    from dppu.torsion.mode import Mode
    from dppu.torsion.nieh_yan import NyVariant
    from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
    print("✓ dppu EC engine modules imported successfully.")
except ImportError as e:
    print(f"✗ Error importing dppu: {e}")
    print("  Ensure script/ (containing dppu/) is on sys.path before importing this module.")
    raise

# =============================================================================
# Constants
# =============================================================================

PI = np.pi
KAPPA = 1.0
L_SCALE = 1.0

_TOPOLOGY_TYPE_MAP = {
    'S3':   TopologyType.S3,
    'T3':   TopologyType.T3,
    'Nil3': TopologyType.NIL3,
}

# Phase classification bounds
R_MIN = 0.01
R_MAX = 1_000_000.0
R_BOUNDARY_THRESHOLD = 0.02

# Ghost detection: small-r probe points (decreasing order).
# Starting from R_BOUNDARY_THRESHOLD down to 0.001 to capture steep divergences.
GHOST_R_SMALL_LIST = (0.02, 0.01, 0.005, 0.002, 0.001)

# α threshold for "α ≈ 0" branch (use Type I/II/III); outside this → new W/G types.
ALPHA_THRESHOLD = 1e-9

# Upper bound for local-minima scan (α ≠ 0 regime).  Using 1e5 rather than R_MAX
# keeps the log-spaced sample dense enough to resolve features near r ~ 1–100.
_COUNT_MINIMA_R_HIGH = 1e5
_COUNT_MINIMA_N      = 300

# =============================================================================
# Engine Registry
# =============================================================================

_MODE_ENUM_MAP = {
    'MX': Mode.MX,
    'VT': Mode.VT,
    'AX': Mode.AX,
}

_NY_ENUM_MAP = {
    'FULL': NyVariant.FULL,
    'TT':   NyVariant.TT,
    'REE':  NyVariant.REE,
}

# =============================================================================
# Potential Function Cache
# =============================================================================
# v4 design: ε and α_W are SymPy symbols inside the engine expression.
# After lambdify, the returned function f(r, V, η, θ_NY, ε, α) can be
# evaluated at any (ε, α) numerically — NO re-run needed.
# Cache key: (topology, mode, ny_variant)  (ε and α NOT in key)

_potential_cache: dict = {}


def get_potential_function(topology: str, mode: str, ny_variant: str,
                           status_callback=None):
    """
    Get (or compute & cache) the potential evaluation function.

    EC design (v4ec): uses UnifiedEngine.  ε and α_W are kept as SymPy
    symbols so a single engine run covers all numerical (ε, α) values.

    EC-specific notes:
      - Nil³: radial symbol is params['R'] (scale param), not params['r'].
      - T³: alpha is not a free symbol (flat → C²≡0); alpha is silently
        ignored in the wrapper.

    Args:
        topology:   'S3', 'T3', or 'Nil3'
        mode:       'MX', 'VT', or 'AX'
        ny_variant: 'FULL', 'TT', or 'REE'
        status_callback: Optional function(str) for progress messages

    Returns:
        Callable: f(r, V_param, eta_param, theta_NY, epsilon, alpha) -> V_eff
    """
    from sympy import lambdify, S

    key = (topology, mode, ny_variant)

    if key in _potential_cache:
        if status_callback:
            status_callback(
                f"Using cached: {topology}×S¹ / {mode} / {ny_variant} "
                f"[ε, α_W evaluated numerically]"
            )
        return _potential_cache[key]

    if status_callback:
        status_callback(
            f"Generating: {topology}×S¹ / {mode} / {ny_variant}…"
            f" (ε/α_W will be numerical after this)"
        )

    from dppu.engine.logger import NullLogger
    cfg = DOFConfig(
        topology=_TOPOLOGY_TYPE_MAP[topology],
        torsion_mode=_MODE_ENUM_MAP[mode],
        ny_variant=_NY_ENUM_MAP[ny_variant],
        enable_squash=True,
        enable_shear=False,
        fiber_mode=FiberMode.NONE,
    )
    engine = UnifiedEngine(cfg, NullLogger())
    engine.run()

    params   = engine.data['params']

    # engine.data['potential'] now uses C²_EC (pipeline step E4.7b → E4.11).
    # No manual correction needed.
    sym_expr = engine.data['potential']

    # Fix L = κ = 1
    sym_expr = sym_expr.subs({params['L']: S.One, params['kappa']: S.One})

    theta_sym   = params['theta_NY']
    V_sym       = params['V']
    eta_sym     = params['eta']
    epsilon_sym = params['epsilon']
    alpha_sym   = params['alpha']

    # Radial symbol: 'r' for S3/T3; 'R' for Nil3 (scale parameter)
    free_syms = sym_expr.free_symbols
    r_sym = params['r'] if params['r'] in free_syms else params['R']

    # Detect which torsion params are free symbols in this mode
    def _is_free(x):
        return hasattr(x, 'is_Symbol') and x.is_Symbol and x in free_syms

    V_is_symbol   = _is_free(V_sym)
    eta_is_symbol = _is_free(eta_sym)
    eps_in_expr   = epsilon_sym in free_syms
    alpha_in_expr = alpha_sym   in free_syms

    # Build lambdify argument list
    lambda_args = [r_sym]
    if V_is_symbol:     lambda_args.append(V_sym)
    if eta_is_symbol:   lambda_args.append(eta_sym)
    lambda_args.append(theta_sym)
    if eps_in_expr:     lambda_args.append(epsilon_sym)
    if alpha_in_expr:   lambda_args.append(alpha_sym)

    raw_func = lambdify(lambda_args, sym_expr, modules='numpy', cse=True)

    # Unified wrapper: f(r, V_param, eta_param, theta_NY, epsilon, alpha)
    def wrapped_func(r, V_param, eta_param, theta_NY, epsilon=0.0, alpha=0.0):
        args = [r]
        if V_is_symbol:   args.append(V_param)
        if eta_is_symbol: args.append(eta_param)
        args.append(theta_NY)
        if eps_in_expr:   args.append(epsilon)
        if alpha_in_expr: args.append(alpha)
        try:
            v = float(np.asarray(raw_func(*args)).flat[0])
            return v if np.isfinite(v) else 1e20
        except Exception:
            return 1e20

    _potential_cache[key] = wrapped_func

    if status_callback:
        status_callback(
            f"✓ Cached: {topology}×S¹ / {mode} / {ny_variant} "
            f"[ε/α_W evaluated numerically at scan time]"
        )

    return wrapped_func


def get_cache_status() -> str:
    """Return current cache contents as a formatted string."""
    if not _potential_cache:
        return "Cache is empty."
    lines = ["Cached configurations  (ε, α_W are evaluated numerically):"]
    for key in sorted(_potential_cache.keys()):
        lines.append(f"  • {key[0]}×S¹ / {key[1]} / {key[2]}")
    return "\n".join(lines)


def clear_cache():
    """Clear the potential function cache."""
    global _potential_cache
    _potential_cache = {}
    print("Cache cleared.")


# =============================================================================
# Weyl Status Helper
# =============================================================================

def get_weyl_c2_status(topology: str, epsilon: float):
    """
    Indicate whether the Weyl scalar C² is zero for the given settings.

    EC physical rules (updated vs v4 LC):
      - T³×S¹:         C²_LC ≡ 0  (flat torus)  → α_W no effect
      - S³×S¹  ε=0:   C²_LC = 0  (conformally flat)  → α_W inactive
      - S³×S¹  ε≠0:   C²_LC ≠ 0  (squashed)  → α_W active
      - Nil³×S¹ any ε: C²_LC = 4/(3R⁴) ≠ 0  (non-conformally flat)
                        → α_W ALWAYS active  [EC paper main result, Theorem 1]

    Returns:
        (is_zero: bool, html_message: str, css_color: str)
    """
    if topology == 'T3':
        msg = (
            "T³ (flat torus): C²_LC ≡ 0 — "
            "<b>α_W has NO effect</b> regardless of ε"
        )
        return True, msg, "#757575"

    if topology == 'Nil3':
        # EC main result: Nil³ is non-conformally flat even at ε = 0.
        # C²_LC(Nil³) = 4/(3R⁴) ≠ 0  (Sec. 5.1, Theorem 1)
        msg = (
            "Nil³ (non-conformally flat): C²_LC = 4/(3R⁴) ≠ 0 — "
            "<b>α_W IS always active</b>"
            "  [EC paper: AX/VT dropout does NOT suppress α at η=V=0]"
        )
        return False, msg, "#1565c0"

    # S³
    if abs(epsilon) < 1e-10:
        msg = (
            "S³, ε = 0 (conformally flat): C²_LC = 0 — "
            "<b>α_W has NO effect</b>"
        )
        return True, msg, "#757575"

    msg = (
        f"S³, ε = {epsilon:.3f} (squashed): C²_LC ≠ 0 — "
        "<b>α_W IS active</b>"
    )
    return False, msg, "#2e7d32"


# =============================================================================
# Phase Classification helpers
# =============================================================================

def _count_local_minima(v_func, r_low=None, r_high=None, n=None) -> int:
    """
    Count local minima of v_func(r) in [r_low, r_high] using log-spaced sampling.

    A point r_k is a local minimum when v(r_k) < v(r_{k-1}) AND v(r_k) < v(r_{k+1}).
    Non-finite values are skipped.

    Used for the α ≠ 0 classification branches (Type G / I-W / II-W).
    """
    if r_low  is None: r_low  = R_BOUNDARY_THRESHOLD
    if r_high is None: r_high = _COUNT_MINIMA_R_HIGH
    if n      is None: n      = _COUNT_MINIMA_N

    r_arr = np.logspace(np.log10(r_low), np.log10(r_high), n)
    vals  = np.empty(n)
    for k, r in enumerate(r_arr):
        try:
            vals[k] = float(np.asarray(v_func(r)).flat[0])
        except Exception:
            vals[k] = np.nan

    count = 0
    for k in range(1, n - 1):
        if (np.isfinite(vals[k - 1]) and np.isfinite(vals[k])
                and np.isfinite(vals[k + 1])
                and vals[k] < vals[k - 1] and vals[k] < vals[k + 1]):
            count += 1
    return count


# =============================================================================
# Phase Classification
# =============================================================================

def classify_phase(potential_func, V_param, eta, theta,
                   epsilon: float = 0.0, alpha: float = 0.0,
                   topology: str = 'S3') -> int:
    """
    Classify stability type of V_eff(r) at given parameters.

    The classification depends on the *effective* sign of α (Weyl coupling).
    When C² ≡ 0 (T³ topology, or ε=0 on S³/Nil³), α·C²=0 regardless of α,
    so the potential is identical to the α=0 case and the standard I/II/III
    logic is used.

    α_eff ≈ 0  (α=0, or C²=0) — standard local classification:
        1: Type I   — local minimum with barrier (V increases toward r=0)
        2: Type II  — rolling (V decreases from r=0 toward minimum)
        3: Type III — no local minimum in physical region

    α_eff > 0  (Ghost regime: V→−∞ as r→0) — minimum-existence check:
        4: Type G     — has local minimum (at best metastable)
        7: Type III-G — no local minimum (unbounded below at r→0)

    α_eff < 0  (Repulsive-wall regime: V→+∞ as r→0) — count local minima:
        6: Type I-W  — two or more local minima
        5: Type II-W — exactly one local minimum
        3: Type III  — no local minimum
    """
    def v(r):
        return potential_func(r, V_param, eta, theta, epsilon, alpha)

    # ── Check whether α is physically active ──────────────────────────────
    # If C² = 0 for this topology+ε, α_W has no effect on V_eff.
    c2_is_zero, _, _ = get_weyl_c2_status(topology, epsilon)
    effective_alpha = 0.0 if c2_is_zero else alpha

    # ── α_eff ≈ 0: existing Type I / II / III logic ───────────────────────
    if abs(effective_alpha) < ALPHA_THRESHOLD:
        try:
            res_min = minimize_scalar(v, bounds=(R_MIN, R_MAX), method='bounded')
        except Exception:
            return 3

        if not res_min.success:
            return 3

        r0 = res_min.x

        if r0 < R_BOUNDARY_THRESHOLD or r0 > R_MAX - R_BOUNDARY_THRESHOLD:
            return 3

        h = 1e-5
        try:
            d2v = (v(r0 + h) - 2 * v(r0) + v(r0 - h)) / h**2
        except Exception:
            return 3

        if d2v <= 0:
            return 3

        # Slope sign near r=0 determines Type I vs II
        r_test = R_MIN * 2
        dr     = R_MIN * 0.1
        try:
            slope_near_zero = (v(r_test + dr) - v(r_test)) / dr
        except Exception:
            slope_near_zero = 0.0

        return 1 if slope_near_zero > 0 else 2

    # ── α_eff > 0: Ghost regime (V→−∞ at r→0) — check for any local minimum
    elif effective_alpha > 0:
        n_min = _count_local_minima(v)
        return 4 if n_min >= 1 else 7  # 7 = Type III-G (unbounded at r→0)

    # ── α_eff < 0: Repulsive-wall regime — count local minima ────────────
    else:
        n_min = _count_local_minima(v)
        if n_min >= 2:
            return 6   # Type I-W  (double well)
        elif n_min == 1:
            return 5   # Type II-W (single well)
        else:
            return 3   # no well


def classify_ghost(potential_func, V_param, eta, theta,
                   epsilon: float = 0.0, alpha: float = 0.0,
                   r_small_list=None, V_cut: float = 1e6) -> tuple:
    """
    Global stability check: does V_eff → −∞ as r → 0?

    This is an independent "global" classification that complements the local
    Type I/II/III classification.  A ghost means any apparent Type I well is
    at best *metastable* — the potential is unbounded below near r=0.

    Algorithm (conservative, purely numerical):
        Evaluate V_eff at a list of decreasing small-r probe points, then
        require ALL three conditions for ghost_flag = True:
          1. Descending tendency: V(r_{k+1}) < V(r_k) for ≥ 70 % of pairs
          2. Sufficiently negative: min(V) < −V_cut
          3. Strong drop (divergent-like): ΔV in the last step > V_cut × 0.1

        Special cases:
          • Any −∞ value at a probe point → ghost immediately (return True)
          • Fewer than 2 finite values     → undetermined (return None)

        Note: r → 0 with V → +∞ (repulsive wall) is NOT a ghost because
        condition 2 (sufficiently negative) will not be satisfied.

    Args:
        potential_func: Callable — same signature as returned by
                        get_potential_function(), i.e.
                        f(r, V_param, eta_param, theta_NY, epsilon, alpha)
        V_param:  Torsion amplitude V
        eta:      Axial torsion η
        theta:    Nieh-Yan coupling θ_NY
        epsilon:  Squashing parameter ε
        alpha:    Weyl coupling α_W
        r_small_list: Probe radii in *decreasing* order.
                      Defaults to module-level GHOST_R_SMALL_LIST.
        V_cut:    Absolute threshold for "sufficiently negative".
                  Overridable via DPPUv2InteractiveViewer.GHOST_V_CUT.

    Returns:
        (ghost_flag: bool | None, ghost_score: float 0.0–1.0)
          True  → ghost confirmed
          None  → undetermined (NaN-dominated; treated as True in display)
          False → no ghost behaviour observed in probe range
    """
    if r_small_list is None:
        r_small_list = GHOST_R_SMALL_LIST

    # ------------------------------------------------------------------
    # Step 1: Evaluate V_eff at each probe point
    # ------------------------------------------------------------------
    vals = []
    for r in r_small_list:
        try:
            v = potential_func(r, V_param, eta, theta, epsilon, alpha)
            v = float(np.asarray(v).flat[0])
        except Exception:
            v = np.nan
        vals.append(v)

    # −∞ at any probe → ghost immediately
    if any(v == -np.inf for v in vals if not np.isnan(v)):
        return True, 1.0

    # Separate finite values for analysis
    finite = [(r, v) for r, v in zip(r_small_list, vals) if np.isfinite(v)]
    if len(finite) < 2:
        return None, 0.0   # undetermined

    v_list = [v for _, v in finite]
    min_v  = min(v_list)

    # ------------------------------------------------------------------
    # Condition 2: Sufficiently negative
    # ------------------------------------------------------------------
    if min_v >= -V_cut:
        return False, 0.0

    # ------------------------------------------------------------------
    # Condition 1: Descending tendency as r decreases
    # (r_small_list is largest→smallest, so v_list[k+1] < v_list[k]
    #  means V keeps falling as r shrinks)
    # ------------------------------------------------------------------
    desc_count = sum(
        1 for i in range(len(v_list) - 1) if v_list[i + 1] < v_list[i]
    )
    desc_ratio = desc_count / (len(v_list) - 1)
    if desc_ratio < 0.7:
        return False, desc_ratio * 0.5

    # ------------------------------------------------------------------
    # Condition 3: Strong (divergent-like) drop in the last step
    # ------------------------------------------------------------------
    _r_prev, v_prev = finite[-2]
    _r_last, v_last = finite[-1]
    delta_v = v_prev - v_last   # positive = still falling as r shrinks
    if delta_v <= V_cut * 0.1:
        return False, desc_ratio * 0.7

    # All three conditions met → ghost confirmed
    score = min(1.0, desc_ratio * 0.5 + min(abs(min_v) / (V_cut * 10), 0.5))
    return True, score


# =============================================================================
# Point Diagnosis  (paper01 slice sense  vs  paper03 full-stationary sense)
# =============================================================================

def diagnose_point(potential_func, r, V_param, eta_param, theta_NY,
                   epsilon=0.0, alpha=0.0, rel_tol=1e-3):
    """
    Diagnose a point (r, V_param, eta_param) in the homogeneous spin-0 sector.

    Computes three numerical partial derivatives via central differences:
      dV/dr, d²V/dr², dV/dη, dV/dV
    and returns a classification in two senses:

    A. Slice extremum (paper01 sense):
       Extremum of V_eff(r) at *fixed* (η, V). Only ∂V/∂r is checked.

    B. Full stationary point (paper03 sense):
       Simultaneous ∂V/∂r = ∂V/∂η = ∂V/∂V = 0.

    Tolerance is scale-adaptive:  tol = rel_tol × max(|V_eff|, 1.0).

    Parameters
    ----------
    potential_func : callable  f(r, V_param, eta_param, theta_NY, epsilon, alpha)
    r             : float      radial coordinate
    V_param       : float      vector torsion amplitude V
    eta_param     : float      axial torsion η
    theta_NY      : float      Nieh-Yan coupling θ
    epsilon       : float      squashing parameter ε  (default 0)
    alpha         : float      Weyl coupling α_W      (default 0)
    rel_tol       : float      relative tolerance for '≈ 0' decisions (default 1e-3)

    Returns
    -------
    dict with keys:
        v_eff          : float  — V_eff at the point
        dv_dr          : float  — ∂V/∂r  (central difference)
        d2v_dr2        : float  — ∂²V/∂r²
        dv_deta        : float  — ∂V/∂η
        dv_dV          : float  — ∂V/∂V
        tol            : float  — tolerance used
        slice_extremum : bool   — |∂V/∂r| < tol
        slice_type     : str    — 'minimum' / 'maximum' / 'flat' / 'not an extremum'
        full_stationary: bool   — all three partials |·| < tol simultaneously
        residual_dir   : str    — 'none' / 'η' / 'V' / 'η and V'
    """
    v0 = float(potential_func(r, V_param, eta_param, theta_NY, epsilon, alpha))

    # Step sizes: relative with floor to avoid zero steps
    hr = max(abs(r),        0.01) * 1e-3
    he = max(abs(eta_param), 0.1) * 1e-3
    hV = max(abs(V_param),   0.1) * 1e-3

    # ∂/∂r  (and ∂²/∂r²)
    try:
        vr_p = float(potential_func(r + hr, V_param,      eta_param,      theta_NY, epsilon, alpha))
        vr_m = float(potential_func(r - hr, V_param,      eta_param,      theta_NY, epsilon, alpha))
        dv_dr   = (vr_p - vr_m) / (2 * hr)
        d2v_dr2 = (vr_p - 2 * v0 + vr_m) / hr ** 2
    except Exception:
        dv_dr = d2v_dr2 = float('nan')

    # ∂/∂η
    try:
        ve_p = float(potential_func(r, V_param, eta_param + he, theta_NY, epsilon, alpha))
        ve_m = float(potential_func(r, V_param, eta_param - he, theta_NY, epsilon, alpha))
        dv_deta = (ve_p - ve_m) / (2 * he)
    except Exception:
        dv_deta = float('nan')

    # ∂/∂V
    try:
        vV_p = float(potential_func(r, V_param + hV, eta_param, theta_NY, epsilon, alpha))
        vV_m = float(potential_func(r, V_param - hV, eta_param, theta_NY, epsilon, alpha))
        dv_dV = (vV_p - vV_m) / (2 * hV)
    except Exception:
        dv_dV = float('nan')

    # Adaptive tolerance
    tol = rel_tol * max(abs(v0), 1.0)

    # A. Slice extremum
    slice_extremum = np.isfinite(dv_dr) and abs(dv_dr) < tol
    if slice_extremum:
        if not np.isfinite(d2v_dr2):
            slice_type = 'inconclusive'
        elif d2v_dr2 > 0:
            slice_type = 'minimum'
        elif d2v_dr2 < 0:
            slice_type = 'maximum'
        else:
            slice_type = 'flat'
    else:
        slice_type = 'not an extremum'

    # B. Full stationary point
    eta_ok = np.isfinite(dv_deta) and abs(dv_deta) < tol
    V_ok   = np.isfinite(dv_dV)   and abs(dv_dV)   < tol
    full_stationary = slice_extremum and eta_ok and V_ok

    # C. Residual slope direction
    if eta_ok and V_ok:
        residual_dir = 'none'
    elif not eta_ok and V_ok:
        residual_dir = 'η'
    elif eta_ok and not V_ok:
        residual_dir = 'V'
    else:
        residual_dir = 'η and V'

    return {
        'v_eff':           v0,
        'dv_dr':           dv_dr,
        'd2v_dr2':         d2v_dr2,
        'dv_deta':         dv_deta,
        'dv_dV':           dv_dV,
        'tol':             tol,
        'slice_extremum':  slice_extremum,
        'slice_type':      slice_type,
        'full_stationary': full_stationary,
        'residual_dir':    residual_dir,
    }


def _format_diagnosis_html(d, r, V_param, eta_param, theta_NY, epsilon, alpha):
    """Format a diagnose_point() result dict as HTML for ipywidgets.HTML."""
    def _fmt(x):
        if not np.isfinite(x):
            return str(x)
        if abs(x) >= 1e4 or (0 < abs(x) < 1e-3):
            return f'{x:.4e}'
        return f'{x:.6g}'

    tol_str = _fmt(d['tol'])

    # Slice extremum colours
    se_color = '#1565c0' if d['slice_extremum'] else '#c62828'
    se_yn    = 'Yes' if d['slice_extremum'] else 'No'
    _st_colors = {
        'minimum':         '#2e7d32',
        'maximum':         '#e65100',
        'flat':            '#757575',
        'inconclusive':    '#757575',
        'not an extremum': '#c62828',
    }
    st_color = _st_colors.get(d['slice_type'], '#555555')

    # Full stationary colours
    fs_color = '#1565c0' if d['full_stationary'] else '#c62828'
    fs_yn    = 'Yes' if d['full_stationary'] else 'No'
    rd_color = '#2e7d32' if d['residual_dir'] == 'none' else '#e65100'
    rd_label = 'none  ✓' if d['residual_dir'] == 'none' else d['residual_dir']

    # Human-readable summary
    if d['full_stationary']:
        msg_en = ('This point <b>IS</b> a full stationary point of the '
                  'homogeneous spin-0 sector.')
    elif d['slice_extremum']:
        msg_en = (
            f'This point is a <b>slice {d["slice_type"]}</b> '
            f'of the fixed-(η,V) radial potential, but <b>NOT</b> a full '
            f'stationary point. Residual slope in: <b>{d["residual_dir"]}</b>.'
        )
    else:
        msg_en = (
            'This point is <b>NOT</b> a slice extremum (∂V/∂r ≠ 0). '
            f'Residual slope in: <b>r and {d["residual_dir"]}</b>.'
        )

    html = (
        '<div style="font-size:12px; font-family:monospace; border:1px solid #ccc;'
        ' padding:8px 12px; border-radius:6px; background:#fafafa; line-height:1.6;">'
        '<b style="font-size:13px;">Point Diagnosis</b>'
        '<hr style="margin:4px 0 6px 0;"/>'
        f'<span style="color:#333;"><b>r</b>={_fmt(r)},  '
        f'<b>V</b>={_fmt(V_param)},  <b>η</b>={_fmt(eta_param)},  '
        f'<b>θ</b>={_fmt(theta_NY)},  <b>ε</b>={_fmt(epsilon)},  '
        f'<b>α</b>={_fmt(alpha)}</span><br/>'
        f'<b>V_eff</b> = {_fmt(d["v_eff"])}'
        f'&emsp;<span style="color:gray; font-size:11px;">[tolerance = {tol_str}]</span>'
        '<hr style="margin:6px 0;"/>'
        '<span style="color:#555;"><b>A. Slice Extremum</b></span>'
        '<span style="font-size:10px; color:gray;"> — fixed-(η,V) radial slice (paper01 sense)</span><br/>'
        f'&nbsp;&nbsp;∂V/∂r&nbsp;&nbsp; = {_fmt(d["dv_dr"])}<br/>'
        f'&nbsp;&nbsp;∂²V/∂r² = {_fmt(d["d2v_dr2"])}<br/>'
        f'&nbsp;&nbsp;<b>Slice extremum?</b> <b style="color:{se_color};">{se_yn}</b>'
        f'&emsp;<b>Type:</b> <b style="color:{st_color};">{d["slice_type"]}</b>'
        '<hr style="margin:6px 0;"/>'
        '<span style="color:#555;"><b>B. Full Stationary Point</b></span>'
        '<span style="font-size:10px; color:gray;"> — homogeneous spin-0 sector (paper03 sense)</span><br/>'
        f'&nbsp;&nbsp;∂V/∂η = {_fmt(d["dv_deta"])}<br/>'
        f'&nbsp;&nbsp;∂V/∂V = {_fmt(d["dv_dV"])}<br/>'
        f'&nbsp;&nbsp;<b>Full stationary?</b> <b style="color:{fs_color};">{fs_yn}</b>'
        '<hr style="margin:6px 0;"/>'
        '<span style="color:#555;"><b>C. Residual Slope Direction:</b></span>'
        f' <b style="color:{rd_color};">{rd_label}</b>'
        '<hr style="margin:6px 0;"/>'
        f'<i style="color:#333;">{msg_en}</i>'
        '</div>'
    )
    return html


# =============================================================================
# Interactive Viewer Class
# =============================================================================

class DPPUv2InteractiveViewer:
    """
    Interactive phase diagram and potential viewer for DPPUv2 (v4 engine).

    New in v4:
      - ε (squashing) slider: deforms S³ and Nil³ geometry
      - α_W (Weyl coupling) slider: enables α·C² term in Lagrangian
      - Weyl status indicator: shows whether C² is zero or active
      - ε and α_W are kept symbolic in engine → ONE run per config

    Usage in notebook::

        # Ensure dppu is on sys.path first:
        import sys, os
        sys.path.insert(0, os.path.abspath('../../'))

        from DPPUv2_interactive_viewer_v4 import DPPUv2InteractiveViewer

        # Optional: customise slider limits
        DPPUv2InteractiveViewer.SLIDER_V_MAX_DEFAULT = 15.0
        DPPUv2InteractiveViewer.SLIDER_ALPHA_W_MAX   = 3.0

        viewer = DPPUv2InteractiveViewer()
        viewer.display()
    """

    # =========================================================================
    # Class-level Configuration (can be modified before instantiation)
    # =========================================================================

    # --- Phase Diagram Range Sliders ---
    SLIDER_V_MAX_MIN     = 1.0
    SLIDER_V_MAX_MAX     = 50.0
    SLIDER_V_MAX_DEFAULT = 10.0

    SLIDER_ETA_MIN_MIN     = -50.0
    SLIDER_ETA_MIN_MAX     = 0.0
    SLIDER_ETA_MIN_DEFAULT = -15.0

    SLIDER_ETA_MAX_MIN     = 0.0
    SLIDER_ETA_MAX_MAX     = 50.0
    SLIDER_ETA_MAX_DEFAULT = 5.0

    # --- Potential Plot Range Sliders (log scale, base-10 exponents) ---
    SLIDER_R_MAX_MIN     = -2    # 10^-2
    SLIDER_R_MAX_MAX     = 6     # 10^6
    SLIDER_R_MAX_DEFAULT = 5.0

    SLIDER_VEFF_MIN_MIN     = 0  # 10^0
    SLIDER_VEFF_MIN_MAX     = 8  # 10^8
    SLIDER_VEFF_MIN_DEFAULT = 1e3

    SLIDER_VEFF_MAX_MIN     = 2  # 10^2
    SLIDER_VEFF_MAX_MAX     = 6  # 10^6
    SLIDER_VEFF_MAX_DEFAULT = 2e3

    # --- Geometric Parameter Sliders (NEW in v4) ---
    # ε: squashing / deformation parameter
    SLIDER_EPSILON_MIN     = -0.9
    SLIDER_EPSILON_MAX     = 3.0
    SLIDER_EPSILON_DEFAULT = 0.0

    # α_W: Weyl coupling  (L = R/2κ² + θ_NY·N + α_W·C²)
    SLIDER_ALPHA_W_MIN     = -5.0
    SLIDER_ALPHA_W_MAX     = 5.0
    SLIDER_ALPHA_W_DEFAULT = 0.0

    # --- Ghost Detection ---
    # Threshold for "sufficiently negative": ghost if min(V_small) < -GHOST_V_CUT
    #
    # Calibration guide (with L_SCALE=κ=1):
    #   GHOST_V_CUT = 100.0  ← default; catches ghost when V_eff < -100 at r~0.001
    #   GHOST_V_CUT = 10.0   ← more sensitive; use if divergence is mild
    #   GHOST_V_CUT = 1e4    ← less sensitive; use if model energy scale is large
    #
    # Use viewer.probe_ghost_values(V, eta) to inspect probe-point values and
    # calibrate this constant for your model.
    GHOST_V_CUT = 100.0

    # =========================================================================
    def __init__(self):
        # --- Core configuration state ---
        self.topology    = 'S3'
        self.torsion_mode = 'MX'
        self.ny_variant  = 'FULL'
        self.theta_NY    = 1.0

        # --- Geometric parameters (v4 new) ---
        self.epsilon     = 0.0   # squashing parameter ε
        self.alpha_weyl  = 0.0   # Weyl coupling α_W

        # --- Phase diagram range ---
        self.V_min   = 0.0
        self.V_max   = 10.0
        self.eta_min = -10.0
        self.eta_max = 10.0

        # --- Potential plot range ---
        self.r_max          = 20.0
        self.Veff_min       = -1e6
        self.Veff_max       = 1e6
        self.r_log_scale    = False
        self.Veff_log_scale = False

        # --- Selected points (up to 3) ---
        self.selected_points = [
            {'V': 4.0, 'eta': -4.0, 'active': True},
            {'V': 4.0, 'eta':  0.0, 'active': True},
            {'V': 4.0, 'eta':  3.0, 'active': True},
        ]
        self.current_point_idx = 0

        # --- Cached phase data ---
        self.phase_data = None
        self.V_grid     = None
        self.eta_grid   = None

        # --- Current potential function ---
        self._current_func   = None
        self._current_config = None

        # --- Matplotlib objects ---
        self.fig         = None
        self.ax_phase    = None
        self.ax_potential = None
        self.click_cid          = None
        self._diagnosis_artists = []   # purple markers added by Diagnose button

        self._create_widgets()

    # =========================================================================
    def _create_widgets(self):
        """Create all ipywidgets controls."""
        style         = {'description_width': '130px'}
        layout_slider = Layout(width='360px')
        layout_drop   = Layout(width='200px')

        # --- Topology / Mode / NY Variant ---
        self.w_topology = widgets.Dropdown(
            options=[('S³×S¹', 'S3'), ('T³×S¹', 'T3'), ('Nil³×S¹', 'Nil3')],
            value='S3',
            description='Topology:',
            style=style, layout=layout_drop
        )

        self.w_torsion_mode = widgets.Dropdown(
            options=[('Mixed (MX)', 'MX'), ('Vector-Trace (VT)', 'VT'), ('Axial (AX)', 'AX')],
            value='MX',
            description='Torsion Mode:',
            style=style, layout=layout_drop
        )

        self.w_ny_variant = widgets.Dropdown(
            options=[('FULL (TT−Ree)', 'FULL'), ('TT only', 'TT'), ('Ree only', 'REE')],
            value='FULL',
            description='NY Variant:',
            style=style, layout=layout_drop
        )

        self.w_theta_NY = widgets.FloatSlider(
            value=1.0, min=0.0, max=10.0, step=0.1,
            description='θ_NY:',
            style=style, layout=layout_slider,
            readout_format='.2f'
        )

        # --- Geometric Parameters (NEW in v4) ---
        self.w_epsilon = widgets.FloatSlider(
            value=self.SLIDER_EPSILON_DEFAULT,
            min=self.SLIDER_EPSILON_MIN,
            max=self.SLIDER_EPSILON_MAX,
            step=0.01,
            description='ε (squashing):',
            style=style, layout=layout_slider,
            readout_format='.3f'
        )

        self.w_alpha_weyl = widgets.FloatSlider(
            value=self.SLIDER_ALPHA_W_DEFAULT,
            min=self.SLIDER_ALPHA_W_MIN,
            max=self.SLIDER_ALPHA_W_MAX,
            step=0.05,
            description='α_W (Weyl):',
            style=style, layout=layout_slider,
            readout_format='.3f'
        )

        # Weyl status indicator — updated dynamically
        _init_c2_zero, _init_msg, _init_color = get_weyl_c2_status('S3', 0.0)
        self.w_weyl_status = widgets.HTML(
            value=(
                f'<span style="color:{_init_color}; font-size:11px;">'
                f'ℹ️ {_init_msg}</span>'
            ),
            layout=Layout(width='100%')
        )

        self.w_geom_box = VBox([
            widgets.HTML(
                '<h5 style="margin:8px 0 2px 0; color:#333;">'
                'Geometric Parameters</h5>'
                '<p style="font-size:10px; color:gray; margin:0 0 4px 0;">'
                'Changes take effect on next "Draw Phase Diagram"</p>'
            ),
            self.w_epsilon,
            self.w_alpha_weyl,
            self.w_weyl_status,
            widgets.HTML(
                '<p style="font-size:10px; color:#555; margin-top:4px;">'
                '  <b>L</b> = R/2κ² + θ_NY·N + <b>α_W·C²</b></p>'
            ),
        ])

        # --- Phase Diagram Range ---
        self.w_V_max = widgets.FloatSlider(
            value=self.SLIDER_V_MAX_DEFAULT,
            min=self.SLIDER_V_MAX_MIN,
            max=self.SLIDER_V_MAX_MAX,
            step=1.0,
            description='V max:',
            style=style, layout=layout_slider
        )

        self.w_eta_min = widgets.FloatSlider(
            value=self.SLIDER_ETA_MIN_DEFAULT,
            min=self.SLIDER_ETA_MIN_MIN,
            max=self.SLIDER_ETA_MIN_MAX,
            step=1.0,
            description='η min:',
            style=style, layout=layout_slider
        )

        self.w_eta_max = widgets.FloatSlider(
            value=self.SLIDER_ETA_MAX_DEFAULT,
            min=self.SLIDER_ETA_MAX_MIN,
            max=self.SLIDER_ETA_MAX_MAX,
            step=1.0,
            description='η max:',
            style=style, layout=layout_slider
        )

        # --- Potential Plot Range ---
        self.w_r_max = widgets.FloatLogSlider(
            value=self.SLIDER_R_MAX_DEFAULT,
            base=10,
            min=self.SLIDER_R_MAX_MIN,
            max=self.SLIDER_R_MAX_MAX,
            step=0.1,
            description='r max:',
            style=style, layout=layout_slider,
            readout_format='.0e'
        )

        self.w_Veff_min = widgets.FloatLogSlider(
            value=self.SLIDER_VEFF_MIN_DEFAULT,
            base=10,
            min=self.SLIDER_VEFF_MIN_MIN,
            max=self.SLIDER_VEFF_MIN_MAX,
            step=0.1,
            description='|V_eff| min:',
            style=style, layout=layout_slider,
            readout_format='.0e'
        )

        self.w_Veff_max = widgets.FloatLogSlider(
            value=self.SLIDER_VEFF_MAX_DEFAULT,
            base=10,
            min=self.SLIDER_VEFF_MAX_MIN,
            max=self.SLIDER_VEFF_MAX_MAX,
            step=0.1,
            description='V_eff max:',
            style=style, layout=layout_slider,
            readout_format='.0e'
        )

        self.w_r_scale = widgets.ToggleButtons(
            options=['Linear', 'Log'],
            value='Linear',
            description='r axis:',
            style=style
        )

        self.w_Veff_scale = widgets.ToggleButtons(
            options=['Linear', 'Symlog'],
            value='Linear',
            description='V_eff axis:',
            style=style
        )

        # --- Point Selection ---
        self.w_point_select = widgets.Dropdown(
            options=[('Point 1 (Black)', 0), ('Point 2 (Green)', 1), ('Point 3 (Red)', 2)],
            value=0,
            description='Active Point:',
            style=style, layout=layout_drop
        )

        self.w_point_V = widgets.FloatSlider(
            value=4.0, min=0.0, max=20.0, step=0.1,
            description='V:',
            style=style, layout=layout_slider
        )

        self.w_point_eta = widgets.FloatSlider(
            value=-4.0, min=-20.0, max=20.0, step=0.1,
            description='η:',
            style=style, layout=layout_slider
        )

        self.w_point_active = widgets.Checkbox(
            value=True,
            description='Show this point',
            style=style
        )

        # --- Status Display ---
        self.w_status = widgets.HTML(
            value='<i style="color:gray;">Ready. Click "Draw Phase Diagram" to start.</i>',
            layout=Layout(width='100%')
        )

        self.w_cache_status = widgets.HTML(
            value=f'<pre style="font-size:11px; color:gray;">{get_cache_status()}</pre>',
            layout=Layout(width='100%')
        )

        # --- Buttons ---
        self.w_draw_button = widgets.Button(
            description='Draw Phase Diagram',
            button_style='primary',
            layout=Layout(width='180px', height='40px')
        )
        self.w_draw_button.on_click(self._on_draw_clicked)

        self.w_update_potential_button = widgets.Button(
            description='Update Potential',
            button_style='success',
            layout=Layout(width='180px', height='40px')
        )
        self.w_update_potential_button.on_click(self._on_update_potential_clicked)

        self.w_clear_cache_button = widgets.Button(
            description='Clear Cache',
            button_style='warning',
            layout=Layout(width='120px', height='30px')
        )
        self.w_clear_cache_button.on_click(self._on_clear_cache_clicked)

        # --- Point Diagnosis (NEW) ---
        self.w_diagnose_r = widgets.BoundedFloatText(
            value=1.0, min=0.001, max=1e6, step=0.001,
            description='r (probe):',
            style=style, layout=Layout(width='260px'),
        )

        self.w_find_extremum_button = widgets.Button(
            description='Find Nearest Extremum',
            button_style='',
            tooltip='Search for slice extremum near current r value',
            layout=Layout(width='195px', height='36px')
        )
        self.w_find_extremum_button.on_click(self._on_find_extremum_clicked)

        self.w_diagnose_button = widgets.Button(
            description='Diagnose This Point',
            button_style='info',
            tooltip='Evaluate gradients at (r, V, η) and classify the point',
            layout=Layout(width='195px', height='36px')
        )
        self.w_diagnose_button.on_click(self._on_diagnose_clicked)

        self.w_diagnosis_output = widgets.HTML(
            value=(
                '<i style="color:gray; font-size:12px;">'
                'Select a point (V, η) above, set r below, '
                'then click "Diagnose This Point".</i>'
            ),
            layout=Layout(width='100%')
        )

        # --- Output Area ---
        self.output = Output()

        # --- Event Handlers ---
        self.w_topology.observe(self._on_topology_or_epsilon_changed, names='value')
        self.w_epsilon.observe(self._on_topology_or_epsilon_changed, names='value')
        self.w_point_select.observe(self._on_point_select_changed, names='value')
        self.w_point_V.observe(self._on_point_slider_changed, names='value')
        self.w_point_eta.observe(self._on_point_slider_changed, names='value')
        self.w_point_active.observe(self._on_point_active_changed, names='value')

    # =========================================================================
    # Status & Cache Helpers
    # =========================================================================

    def _set_status(self, msg: str, color: str = 'black'):
        self.w_status.value = f'<span style="color:{color};">{msg}</span>'

    def _update_cache_display(self):
        self.w_cache_status.value = (
            f'<pre style="font-size:11px; color:gray;">{get_cache_status()}</pre>'
        )

    def _update_weyl_status(self):
        """Refresh the Weyl C² status indicator based on current widget values."""
        topo = self.w_topology.value
        eps  = self.w_epsilon.value
        _, msg, color = get_weyl_c2_status(topo, eps)
        self.w_weyl_status.value = (
            f'<span style="color:{color}; font-size:11px;">ℹ️ {msg}</span>'
        )

    # =========================================================================
    # Ghost Probe Diagnostic
    # =========================================================================

    def probe_ghost_values(self, V_param=None, eta=None,
                           theta=None, epsilon=None, alpha=None,
                           r_list=None):
        """
        Diagnostic: print V_eff at each small-r probe point and show which
        classify_ghost() conditions are met.

        Useful for calibrating GHOST_V_CUT when ghost detection is not firing.

        Args:
            V_param: V torsion amplitude (default: current Point 1)
            eta:     Axial torsion η     (default: current Point 1)
            theta:   θ_NY               (default: current widget value)
            epsilon: ε                  (default: current widget value)
            alpha:   α_W                (default: current widget value)
            r_list:  Custom probe list  (default: GHOST_R_SMALL_LIST)

        Usage in notebook::

            viewer.probe_ghost_values(V_param=3.0, eta=-4.0)
            # or: probe Point 1 with current widget settings
            viewer.probe_ghost_values()
        """
        # Defaults from current state
        if V_param is None:
            V_param = self.selected_points[0]['V']
        if eta is None:
            eta = self.selected_points[0]['eta']
        if theta is None:
            theta = self.w_theta_NY.value
        if epsilon is None:
            epsilon = self.w_epsilon.value
        if alpha is None:
            alpha = self.w_alpha_weyl.value
        if r_list is None:
            r_list = GHOST_R_SMALL_LIST

        func = self._get_potential_func()
        V_cut = self.GHOST_V_CUT

        print("=" * 58)
        print(f"Ghost Probe: V={V_param:.3f}, η={eta:.3f}, θ={theta:.3f}")
        print(f"             ε={epsilon:.3f}, α_W={alpha:.3f}")
        print(f"  GHOST_V_CUT = {V_cut}")
        print("-" * 58)
        print(f"  {'r':>8}   {'V_eff':>14}   notes")
        print("-" * 58)

        vals = []
        for r in r_list:
            try:
                v = func(r, V_param, eta, theta, epsilon, alpha)
                v = float(np.asarray(v).flat[0])
            except Exception as ex:
                v = np.nan
                print(f"  {r:>8.4f}   {'NaN':>14}   (exception: {ex})")
                vals.append(np.nan)
                continue
            vals.append(v)
            note = ""
            if np.isneginf(v):
                note = "← −∞ detected!"
            elif not np.isfinite(v):
                note = "← non-finite"
            elif v < -V_cut:
                note = f"← below −V_cut ({-V_cut})"
            print(f"  {r:>8.4f}   {v:>14.4g}   {note}")

        print("-" * 58)

        # Run full classify_ghost for diagnosis
        ghost_flag, score = classify_ghost(
            func, V_param, eta, theta,
            epsilon=epsilon, alpha=alpha, V_cut=V_cut
        )

        finite = [(r, v) for r, v in zip(r_list, vals) if np.isfinite(v)]
        if len(finite) >= 2:
            v_list = [v for _, v in finite]
            min_v = min(v_list)
            desc_count = sum(
                1 for i in range(len(v_list)-1) if v_list[i+1] < v_list[i]
            )
            desc_ratio = desc_count / (len(v_list) - 1)
            delta_v = v_list[-2] - v_list[-1] if len(v_list) >= 2 else 0.0

            cond1 = "✓" if desc_ratio >= 0.7 else "✗"
            cond2 = "✓" if min_v < -V_cut else "✗"
            cond3 = "✓" if delta_v > V_cut * 0.1 else "✗"
            print(f"  Cond 1 (desc≥70%): {cond1}  desc_ratio={desc_ratio:.2f}")
            print(f"  Cond 2 (min<-V_cut): {cond2}  min={min_v:.4g}  threshold={-V_cut}")
            print(f"  Cond 3 (Δv>{V_cut*0.1:.3g}): {cond3}  Δv={delta_v:.4g}  (last step)")
        print("-" * 58)
        flag_str = {True: "GHOST", False: "no ghost", None: "undetermined"}
        print(f"  → ghost_flag = {flag_str.get(ghost_flag, str(ghost_flag))},  score = {score:.3f}")
        if ghost_flag is not True:
            print(f"  Tip: if ghost is visible in the plot, try lowering GHOST_V_CUT")
            print(f"       DPPUv2InteractiveViewer.GHOST_V_CUT = <value smaller than {abs(min_v):.3g}>")
        print("=" * 58)

    # =========================================================================
    # Widget Event Handlers
    # =========================================================================

    def _on_topology_or_epsilon_changed(self, change):
        """Update Weyl status indicator when topology or ε changes."""
        self._update_weyl_status()

    def _on_point_select_changed(self, change):
        idx = change['new']
        self.current_point_idx = idx
        pt = self.selected_points[idx]

        # Update sliders without cascading observer callbacks
        for w in (self.w_point_V, self.w_point_eta):
            w.unobserve(self._on_point_slider_changed, names='value')
        self.w_point_active.unobserve(self._on_point_active_changed, names='value')

        self.w_point_V.value     = pt['V']
        self.w_point_eta.value   = pt['eta']
        self.w_point_active.value = pt['active']

        for w in (self.w_point_V, self.w_point_eta):
            w.observe(self._on_point_slider_changed, names='value')
        self.w_point_active.observe(self._on_point_active_changed, names='value')

    def _on_point_slider_changed(self, change):
        idx = self.current_point_idx
        self.selected_points[idx]['V']   = self.w_point_V.value
        self.selected_points[idx]['eta'] = self.w_point_eta.value
        self._update_point_markers()

    def _on_point_active_changed(self, change):
        idx = self.current_point_idx
        self.selected_points[idx]['active'] = change['new']
        self._update_point_markers()
        self._update_potential_plot()

    def _on_draw_clicked(self, b):
        self._read_widget_values()
        self._compute_phase_diagram()
        self._draw_all()
        self._update_cache_display()

    def _on_update_potential_clicked(self, b):
        self._read_widget_values()
        self._update_potential_plot()

    def _on_clear_cache_clicked(self, b):
        clear_cache()
        self._current_func   = None
        self._current_config = None
        self._update_cache_display()
        self._set_status("Cache cleared.", color='orange')

    def _on_find_extremum_clicked(self, b):
        """
        Search for the nearest slice extremum (local minimum) around the
        current w_diagnose_r value for the active point.
        """
        func  = self._get_potential_func()
        idx   = self.current_point_idx
        pt    = self.selected_points[idx]
        theta = self.w_theta_NY.value
        eps   = self.w_epsilon.value
        alp   = self.w_alpha_weyl.value

        r_cur = self.w_diagnose_r.value
        r_lo  = max(R_MIN, r_cur * 0.05)
        r_hi  = min(max(self.r_max, 100.0), r_cur * 20.0)
        if r_hi <= r_lo:
            r_lo, r_hi = R_MIN, max(self.r_max, 100.0)

        def v1d(r_):
            return func(r_, pt['V'], pt['eta'], theta, eps, alp)

        try:
            res = minimize_scalar(
                v1d, bounds=(r_lo, r_hi), method='bounded',
                options={'xatol': 1e-6}
            )
            if not res.success:
                self._set_status(
                    'Could not find slice extremum in the search range.',
                    color='orange'
                )
                return

            r_found = float(res.x)
            h       = r_found * 1e-3
            v0      = float(v1d(r_found))
            d2v     = (float(v1d(r_found + h)) - 2 * v0
                       + float(v1d(r_found - h))) / h ** 2

            if r_lo * 1.01 < r_found < r_hi * 0.99 and d2v > 0:
                # FloatLogSlider.value holds the actual float value
                self.w_diagnose_r.value = r_found
                self._set_status(
                    f'Slice extremum found: r = {r_found:.5g}  '
                    f'(d\u00b2V/dr\u00b2 = {d2v:.4g})',
                    color='green'
                )
                # 判定ロジックを走らせて紫点を再描画
                self._on_diagnose_clicked(None)
            else:
                self._set_status(
                    f'Potential converged to boundary r = {r_found:.5g} '
                    '\u2014 no interior extremum in range.',
                    color='orange'
                )
        except Exception as ex:
            self._set_status(f'Find extremum error: {ex}', color='red')

    def _on_diagnose_clicked(self, b):
        """
        Run point diagnosis for (r, V, η) of the active point and display
        the result in the diagnosis panel. Also marks the point on the
        potential plot with a purple marker.
        """
        func  = self._get_potential_func()
        idx   = self.current_point_idx
        pt    = self.selected_points[idx]
        r     = self.w_diagnose_r.value
        theta = self.w_theta_NY.value
        eps   = self.w_epsilon.value
        alp   = self.w_alpha_weyl.value

        try:
            result   = diagnose_point(
                func, r, pt['V'], pt['eta'], theta,
                epsilon=eps, alpha=alp
            )
            html_str = _format_diagnosis_html(
                result, r, pt['V'], pt['eta'], theta, eps, alp
            )
            self.w_diagnosis_output.value = html_str

            # Mark the diagnosis point on the potential axe
            if self.ax_potential is not None and self.fig is not None:
                for artist in self._diagnosis_artists:
                    try:
                        artist.remove()
                    except Exception:
                        pass
                self._diagnosis_artists = []

                v_at_r = result['v_eff']
                ylo, yhi = self.ax_potential.get_ylim()
                v_display = float(np.clip(v_at_r, ylo, yhi))

                vline = self.ax_potential.axvline(
                    r, color='purple', linestyle=':', linewidth=1.5, alpha=0.7
                )
                pt_marker, = self.ax_potential.plot(
                    [r], [v_display], 'D',
                    color='purple', markersize=8, zorder=15,
                    label=f'Diagnosis r={r:.4g}'
                )
                self._diagnosis_artists = [vline, pt_marker]
                self.ax_potential.legend(loc='best', fontsize=9)
                self.fig.canvas.draw_idle()

        except Exception as ex:
            self.w_diagnosis_output.value = (
                f'<span style="color:red;">Diagnosis error: {ex}</span>'
            )

    # =========================================================================
    # State Management
    # =========================================================================

    def _read_widget_values(self):
        """Synchronise internal state from current widget values."""
        self.topology     = self.w_topology.value
        self.torsion_mode = self.w_torsion_mode.value
        self.ny_variant   = self.w_ny_variant.value
        self.theta_NY     = self.w_theta_NY.value

        # Geometric parameters (v4)
        self.epsilon     = self.w_epsilon.value
        self.alpha_weyl  = self.w_alpha_weyl.value

        self.V_min   = 0.0
        self.V_max   = self.w_V_max.value
        self.eta_min = self.w_eta_min.value
        self.eta_max = self.w_eta_max.value

        self.r_max          = self.w_r_max.value
        self.Veff_min       = -self.w_Veff_min.value
        self.Veff_max       = self.w_Veff_max.value
        self.r_log_scale    = (self.w_r_scale.value == 'Log')
        self.Veff_log_scale = (self.w_Veff_scale.value == 'Symlog')

    def _get_potential_func(self):
        """Get the potential function (with caching)."""
        config = (self.topology, self.torsion_mode, self.ny_variant)
        if self._current_config != config:
            self._current_func = get_potential_function(
                self.topology,
                self.torsion_mode,
                self.ny_variant,
                status_callback=self._set_status
            )
            self._current_config = config
        return self._current_func

    # =========================================================================
    # Phase Diagram Computation
    # =========================================================================

    def _compute_phase_diagram(self):
        """Compute 100×100 phase classification grid."""
        n_points = 100

        V_range   = np.linspace(self.V_min, self.V_max, n_points)
        eta_range = np.linspace(self.eta_min, self.eta_max, n_points)

        self.V_grid, self.eta_grid = np.meshgrid(V_range, eta_range)
        self.phase_data = np.zeros_like(self.V_grid)

        func = self._get_potential_func()

        # Build title suffix for status messages
        status_str = (
            f"{self.topology}×S¹ / {self.torsion_mode} / {self.ny_variant} / "
            f"ε={self.epsilon:.3f} / α_W={self.alpha_weyl:.3f}"
        )
        self._set_status(f"Computing phase diagram for {status_str}…", color='blue')

        for i in range(n_points):
            for j in range(n_points):
                self.phase_data[i, j] = classify_phase(
                    func,
                    self.V_grid[i, j],
                    self.eta_grid[i, j],
                    self.theta_NY,
                    epsilon=self.epsilon,
                    alpha=self.alpha_weyl,
                    topology=self.topology,
                )

            if (i + 1) % 10 == 0:
                pct = 100 * (i + 1) / n_points
                self._set_status(
                    f"Computing phase diagram… {pct:.0f}%", color='blue'
                )

        self._set_status(
            f"✓ Phase diagram computed: {status_str} / θ_NY={self.theta_NY:.2f}",
            color='green'
        )

    # =========================================================================
    # Drawing
    # =========================================================================

    def _draw_all(self):
        """Draw both the phase diagram and potential plot."""
        with self.output:
            clear_output(wait=True)

            if self.fig is not None:
                plt.close(self.fig)

            self.fig, (self.ax_phase, self.ax_potential) = plt.subplots(
                1, 2, figsize=(14, 5), constrained_layout=True
            )

            self._draw_phase_diagram()
            self._draw_potential_plot()

            if self.click_cid is not None:
                self.fig.canvas.mpl_disconnect(self.click_cid)
            self.click_cid = self.fig.canvas.mpl_connect(
                'button_press_event', self._on_click
            )

            plt.show()

    def _draw_phase_diagram(self):
        ax = self.ax_phase
        ax.clear()

        if self.phase_data is None:
            ax.text(0.5, 0.5, 'Click "Draw Phase Diagram" to compute',
                    ha='center', va='center', transform=ax.transAxes)
            return

        # ── Colours ──────────────────────────────────────────────────────────
        # α ≈ 0  (standard)
        color_I    = '#90EE90'   # light green  — Type I    stable well
        color_II   = '#FFFF99'   # light yellow — Type II   rolling
        color_III  = '#E8E8E8'   # light gray   — Type III  no well (any α)
        # α > 0  (ghost / unbounded below)
        color_G    = '#FFCCCC'   # light red    — Type G    metastable well
        color_IIIG = '#E8E8E8'   # light gray   — Type III-G no well + unbounded (same as III, hatched)
        # α < 0  (repulsive wall)
        color_IIW  = '#FFE0B2'   # light orange — Type II-W single well
        color_IW   = '#BBDEFB'   # light blue   — Type I-W  double well
        color_err  = '#FFFFFF'   # white        — error / unclassified

        # Integer → colour mapping
        #  0: err  1: I  2: II  3: III  4: G  5: II-W  6: I-W  7: III-G
        cmap   = mcolors.ListedColormap(
            [color_err, color_I, color_II, color_III, color_G, color_IIW, color_IW, color_IIIG]
        )
        bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]
        norm   = mcolors.BoundaryNorm(bounds, cmap.N)

        ax.pcolormesh(self.V_grid, self.eta_grid, self.phase_data,
                      cmap=cmap, norm=norm, shading='auto')

        # ── Type III-G hatching overlay (α>0, no local minimum) ─────────────
        if np.any(self.phase_data == 7):
            _iiig_bin = (self.phase_data == 7).astype(float)
            _n_before = len(ax.collections)
            ax.contourf(self.V_grid, self.eta_grid, _iiig_bin,
                        levels=[0.5, 1.5], colors=['none'], hatches=['///'])
            # ax.collections is stable across matplotlib versions (unlike cf.collections)
            for _coll in ax.collections[_n_before:]:
                _coll.set_edgecolor('#808080')
                _coll.set_linewidth(0.3)

        # ── Phase boundary contour lines ────────────────────────────────────
        # Use effective α: if C²=0 (T3 or ε=0), α has no physical effect.
        _c2_zero, _, _ = get_weyl_c2_status(self.topology, self.epsilon)
        alpha_now = 0.0 if _c2_zero else self.alpha_weyl
        if abs(alpha_now) < ALPHA_THRESHOLD:
            # Standard: I/II at 1.5, II/III at 2.5
            ax.contour(self.V_grid, self.eta_grid, self.phase_data,
                       levels=[1.5, 2.5], colors=['darkgreen', 'darkorange'],
                       linewidths=1.5)
        elif alpha_now > 0:
            # Ghost: G(4) / III-G(7) boundary at 5.5 (only G and III-G appear)
            ax.contour(self.V_grid, self.eta_grid, self.phase_data,
                       levels=[5.5], colors=['darkred'],
                       linewidths=1.5)
        else:
            # Repulsive wall: III(3)/II-W(5) at 4.0, II-W(5)/I-W(6) at 5.5
            ax.contour(self.V_grid, self.eta_grid, self.phase_data,
                       levels=[4.0, 5.5], colors=['darkorange', 'darkblue'],
                       linewidths=1.5)

        self._draw_point_markers(ax)

        ax.set_xlabel(r'$V$ (Vector Torsion)', fontsize=11)
        ax.set_ylabel(r'$\eta$ (Axial Torsion)', fontsize=11)

        # Title: show all relevant parameters
        c2_note = ""
        is_zero, _, _ = get_weyl_c2_status(self.topology, self.epsilon)
        if not is_zero:
            c2_note = f", α_W={self.alpha_weyl:.3f}·C²"
        title = (
            f"{self.topology}×S¹ Phase Diagram  "
            f"(ε={self.epsilon:.3f}{c2_note})\n"
            f"{self.torsion_mode} / {self.ny_variant},  θ_NY={self.theta_NY:.2f}"
        )
        ax.set_title(title, fontsize=11)
        ax.set_xlim(self.V_min, self.V_max)
        ax.set_ylim(self.eta_min, self.eta_max)
        ax.grid(True, alpha=0.3, linestyle='--')

        # ── Legend: show only the types relevant to current α ────────────────
        if abs(alpha_now) < ALPHA_THRESHOLD:
            legend_elements = [
                Patch(facecolor=color_I,   edgecolor='darkgreen',  label='Type I (Stable Well)'),
                Patch(facecolor=color_II,  edgecolor='darkorange', label='Type II (Rolling)'),
                Patch(facecolor=color_III, edgecolor='gray',       label='Type III (No Well)'),
            ]
        elif alpha_now > 0:
            legend_elements = [
                Patch(facecolor=color_G,    edgecolor='darkred',    label='Type G (Metastable, V→−∞ at r→0)'),
                Patch(facecolor=color_IIIG, edgecolor='gray',       label='Type III-G (No Well, V→−∞ at r→0)',
                      hatch='///'),
            ]
        else:
            legend_elements = [
                Patch(facecolor=color_IW,  edgecolor='darkblue',   label='Type I-W (Double Well, α<0)'),
                Patch(facecolor=color_IIW, edgecolor='darkorange', label='Type II-W (Single Well, α<0)'),
                Patch(facecolor=color_III, edgecolor='gray',       label='Type III (No Well)'),
            ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

    def _draw_point_markers(self, ax):
        markers = ['o', 's', '^']
        colors  = ['black', 'green', 'red']
        sizes   = [100, 80, 80]
        for i, (pt, m, c, s) in enumerate(
                zip(self.selected_points, markers, colors, sizes)):
            if pt['active']:
                ax.scatter(pt['V'], pt['eta'],
                           marker=m, s=s, facecolors='none',
                           edgecolors=c, linewidths=2, zorder=10)

    def _update_point_markers(self):
        if self.ax_phase is not None and self.phase_data is not None:
            self._draw_phase_diagram()
            self.fig.canvas.draw_idle()

    def _draw_potential_plot(self):
        ax = self.ax_potential
        ax.clear()
        self._diagnosis_artists = []   # cleared by ax.clear(); reset reference list

        func = self._get_potential_func()

        if self.r_log_scale:
            r_arr = np.logspace(np.log10(0.01), np.log10(self.r_max), 1000)
        else:
            r_arr = np.linspace(0.01, self.r_max, 1000)

        colors     = ['black', 'green', 'red']
        linestyles = ['-', '--', ':']
        labels     = ['Point 1', 'Point 2', 'Point 3']

        has_plot = False
        for i, pt in enumerate(self.selected_points):
            if pt['active']:
                V_arr = np.array([
                    func(r, pt['V'], pt['eta'], self.theta_NY,
                         self.epsilon, self.alpha_weyl)
                    for r in r_arr
                ])
                phase_type = classify_phase(
                    func, pt['V'], pt['eta'], self.theta_NY,
                    self.epsilon, self.alpha_weyl,
                    topology=self.topology,
                )
                _TYPE_LABELS = {
                    0: 'Err',
                    1: 'Type I',
                    2: 'Type II',
                    3: 'Type III',
                    4: 'Type G',
                    5: 'Type II-W',
                    6: 'Type I-W',
                    7: 'Type III-G',
                }
                type_str = _TYPE_LABELS.get(phase_type, f'Type {phase_type}')

                ax.plot(r_arr, V_arr,
                        color=colors[i], linestyle=linestyles[i],
                        linewidth=2,
                        label=f"{labels[i]}: V={pt['V']:.1f}, η={pt['eta']:.1f} ({type_str})")
                has_plot = True

        if not has_plot:
            ax.text(0.5, 0.5, 'No points selected',
                    ha='center', va='center', transform=ax.transAxes)
            return

        ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)

        if self.r_log_scale:
            ax.set_xscale('log')
        if self.Veff_log_scale:
            ax.set_yscale('symlog', linthresh=100)

        ax.set_xlim(0.01 if self.r_log_scale else 0, self.r_max)
        ax.set_ylim(self.Veff_min, self.Veff_max)

        ax.set_xlabel(r'$r$ (Radius)', fontsize=11)
        ax.set_ylabel(r'$V_{\mathrm{eff}}(r)$', fontsize=11)

        # Potential title: show ε and α_W if relevant
        is_zero, _, _ = get_weyl_c2_status(self.topology, self.epsilon)
        if is_zero:
            geom_note = f"ε={self.epsilon:.3f}  (α_W inactive)"
        else:
            geom_note = f"ε={self.epsilon:.3f}, α_W={self.alpha_weyl:.3f}"
        ax.set_title(
            f"Effective Potential  θ_NY={self.theta_NY:.2f}\n({geom_note})",
            fontsize=11
        )

        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3, linestyle='--')

    def _update_potential_plot(self):
        if self.ax_potential is not None:
            self._read_widget_values()
            self._draw_potential_plot()
            self.fig.canvas.draw_idle()

    # =========================================================================
    # Click Handler
    # =========================================================================

    def _on_click(self, event):
        if event.button != 1:
            return

        # ── Click on phase diagram: select (V, η) point ─────────────────────
        if event.inaxes == self.ax_phase:
            V_clicked   = event.xdata
            eta_clicked = event.ydata
            if V_clicked is None or eta_clicked is None:
                return

            V_clicked   = np.clip(V_clicked,   self.V_min,   self.V_max)
            eta_clicked = np.clip(eta_clicked, self.eta_min, self.eta_max)

            idx = self.current_point_idx
            self.selected_points[idx]['V']      = V_clicked
            self.selected_points[idx]['eta']    = eta_clicked
            self.selected_points[idx]['active'] = True

            for w in (self.w_point_V, self.w_point_eta):
                w.unobserve(self._on_point_slider_changed, names='value')
            self.w_point_active.unobserve(self._on_point_active_changed, names='value')

            self.w_point_V.value      = V_clicked
            self.w_point_eta.value    = eta_clicked
            self.w_point_active.value = True

            for w in (self.w_point_V, self.w_point_eta):
                w.observe(self._on_point_slider_changed, names='value')
            self.w_point_active.observe(self._on_point_active_changed, names='value')

            self._draw_phase_diagram()
            self._draw_potential_plot()
            self.fig.canvas.draw_idle()

        # ── Click on potential plot: set diagnosis r ──────────────────────────
        elif event.inaxes == self.ax_potential:
            r_clicked = event.xdata
            if r_clicked is None or r_clicked <= 0:
                return
            # Clamp to BoundedFloatText limits
            r_clicked = float(np.clip(r_clicked, 0.001, 1e6))
            self.w_diagnose_r.value = r_clicked
            # Run diagnosis immediately so the purple marker appears on click
            self._on_diagnose_clicked(None)

    # =========================================================================
    # Display
    # =========================================================================

    def display(self):
        """
        Render the interactive viewer in the Jupyter notebook.

        Call this method after (optionally) adjusting class-level defaults::

            DPPUv2InteractiveViewer.SLIDER_EPSILON_MAX  = 2.0
            DPPUv2InteractiveViewer.SLIDER_ALPHA_W_MAX  = 3.0
            viewer = DPPUv2InteractiveViewer()
            viewer.display()
        """
        from IPython.display import clear_output
        clear_output(wait=False)

        # ── Top-left: Configuration ────────────────────────────────────────
        config_box = VBox([
            widgets.HTML('<h4 style="margin-bottom:4px;">Configuration</h4>'),
            HBox([self.w_topology, self.w_torsion_mode, self.w_ny_variant]),
            self.w_theta_NY,
            self.w_geom_box,
        ])

        # ── Top-right: Status ─────────────────────────────────────────────
        status_box = VBox([
            widgets.HTML('<h4 style="margin-bottom:4px;">Status</h4>'),
            self.w_status,
            HBox([self.w_clear_cache_button]),
            self.w_cache_status,
        ])

        top_section = HBox([config_box, status_box])

        # ── Bottom-left: Phase Diagram Controls ────────────────────────────
        phase_controls = VBox([
            widgets.HTML(
                '<h4 style="border-bottom:2px solid #4CAF50; padding-bottom:5px;">'
                'Phase Diagram Controls</h4>'
            ),
            widgets.HTML('<h5>Range</h5>'),
            self.w_V_max,
            HBox([self.w_eta_min, self.w_eta_max]),
            widgets.HTML('<h5 style="margin-top:10px;">Point Selection</h5>'),
            widgets.HTML(
                '<p style="font-size:11px; color:gray;">'
                'Click on phase diagram or use sliders below</p>'
            ),
            self.w_point_select,
            HBox([self.w_point_V, self.w_point_eta]),
            self.w_point_active,
            widgets.HTML('<div style="margin-top:15px;"></div>'),
            self.w_draw_button,
        ])

        # ── Bottom-right: Potential Plot Controls ─────────────────────────
        potential_controls = VBox([
            widgets.HTML(
                '<h4 style="border-bottom:2px solid #2196F3; padding-bottom:5px;">'
                'Potential Plot Controls</h4>'
            ),
            widgets.HTML('<h5>Range</h5>'),
            self.w_r_max,
            HBox([self.w_Veff_min, self.w_Veff_max]),
            widgets.HTML('<h5 style="margin-top:10px;">Axis Scale</h5>'),
            HBox([self.w_r_scale, self.w_Veff_scale]),
            widgets.HTML('<div style="margin-top:15px;"></div>'),
            self.w_update_potential_button,
        ])

        bottom_section = HBox([phase_controls, potential_controls])

        # ── Point Diagnosis Controls (above graph) ───────────────────────
        diagnosis_controls = VBox([
            widgets.HTML(
                '<h4 style="border-bottom:2px solid #9C27B0; padding-bottom:5px;">'
                'Point Diagnosis</h4>'
                '<p style="font-size:11px; color:gray; margin:2px 0 6px 0;">'
                'Diagnose <b>Slice Extremum</b> (fixed-(η,V) sense — paper01) '
                'vs <b>Full Stationary Point</b> (spin-0 sector sense — paper03).</p>'
                '<p style="font-size:11px; color:#555; margin:0 0 6px 0;">'
                '▶ Click on the <b>potential plot</b> to set <i>r</i> and run diagnosis instantly.  '
                'Or type a value in the box and click <b>Diagnose This Point</b>.<br/>'
                '⚠️ Diagnosis applies to the <b>Active Point only</b> '
                '(currently selected via the “Active Point” dropdown above). '
                'Points 2 and 3 are not diagnosed.</p>'
            ),
            HBox([
                self.w_diagnose_r,
                self.w_find_extremum_button,
                self.w_diagnose_button,
            ]),
        ], layout=Layout(width='100%', padding='4px 0'))

        controls = VBox([top_section, bottom_section, diagnosis_controls])
        display(controls)
        display(self.output)

        # Initial placeholder plot
        # NOTE: display(self.w_diagnosis_output) must come AFTER this block.
        # In %matplotlib widget mode, calling display() between display(self.output)
        # and the 'with self.output: plt.show()' block causes the figure canvas
        # to be pushed to both the Output widget and the cell output (double render).
        with self.output:
            self.fig, (self.ax_phase, self.ax_potential) = plt.subplots(
                1, 2, figsize=(14, 5), constrained_layout=True
            )
            self.ax_phase.text(
                0.5, 0.5, 'Click "Draw Phase Diagram" to start',
                ha='center', va='center',
                transform=self.ax_phase.transAxes, fontsize=12, color='gray'
            )
            self.ax_phase.set_title('Phase Diagram')
            self.ax_potential.text(
                0.5, 0.5, 'Potential will appear here',
                ha='center', va='center',
                transform=self.ax_potential.transAxes, fontsize=12, color='gray'
            )
            self.ax_potential.set_title('Effective Potential')

            self.click_cid = self.fig.canvas.mpl_connect(
                'button_press_event', self._on_click
            )
            plt.show()

        # ── Diagnosis result panel (below graph) ────────────────────────
        # Placed after the 'with self.output:' block to avoid double-render.
        display(self.w_diagnosis_output)


# =============================================================================
# Usage Example
# =============================================================================
# In a Jupyter notebook:
#
#   import sys, os
#   sys.path.insert(0, os.path.abspath('../../'))   # add _DPPUv2_Phase2/
#
#   from DPPUv2_interactive_viewer_v4 import DPPUv2InteractiveViewer, clear_cache
#
#   # Optional: customise slider defaults
#   DPPUv2InteractiveViewer.SLIDER_V_MAX_DEFAULT    = 15.0
#   DPPUv2InteractiveViewer.SLIDER_ETA_MIN_DEFAULT  = -20.0
#   DPPUv2InteractiveViewer.SLIDER_EPSILON_MAX      = 2.0
#   DPPUv2InteractiveViewer.SLIDER_ALPHA_W_MIN      = -3.0
#   DPPUv2InteractiveViewer.SLIDER_ALPHA_W_MAX      = 3.0
#
#   viewer = DPPUv2InteractiveViewer()
#   viewer.display()
