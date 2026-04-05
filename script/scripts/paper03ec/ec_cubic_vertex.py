"""
EC Cubic Vertex Analysis
=========================

Analyzes the cubic vertices (third-order derivatives) of the EC-Weyl
effective potential V_eff_EC across topologies.

Key results:
  1. theta-cubic (NY channel): preserved from LC in EC
     d3V/d(theta)d(eta)d(V) is unchanged by the EC coupling

  2. alpha-cubic (new EC channel):
     d3V/d(alpha)d(V)d(eta) != 0
     This is the EC-specific homogeneous EFT cubic coefficient

  3. Topology comparison:
     S3, T3, Nil3 all exhibit the same theta-cubic preservation
     while the detailed mixed cubic coefficients differ by topology

  4. Pontryagin coefficient:
     the MIXING-mode cubic coefficient is alpha-independent in the T3 check

Author: Muacca
Date: 2026-03-30
"""

import sys
import os

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_LIB_DIR    = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", ".."))
_DATA_DIR   = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "..", "data"))
sys.path.insert(0, _LIB_DIR)
os.makedirs(_DATA_DIR, exist_ok=True)

from dppu.utils.tee_logger import setup_log, teardown_log
setup_log(__file__, log_dir=_DATA_DIR)

from sympy import cancel, diff, factor, S

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.action.ec_action import build_veff_ec

print("=" * 70)
print("EC Cubic Vertex Analysis")
print("=" * 70)
print()

# ── Helper: V_eff for each topology ──────────────────────────────────────────

def get_veff(topo_type, name):
    print(f"  [Building {name} MX engine...]", flush=True)
    veff_ec, veff_lc, params = build_veff_ec(
        topo_type,
        torsion_mode=Mode.MX,
        ny_variant=NyVariant.FULL,
        L_val=1,
        kappa_val=1,
    )
    print("  Done")
    return veff_ec, veff_lc, params


# ── Part 1: S3 cubic vertex analysis ─────────────────────────────────────────
print("-" * 70)
print("Part 1: S3 cubic vertex analysis")
print("-" * 70)
print()

veff_ec_s3, veff_lc_s3, params_s3 = get_veff(TopologyType.S3, "S3")

r_s3     = params_s3['r']
eta_s3   = params_s3['eta']
V_s3     = params_s3['V']
alpha_s3 = params_s3['alpha']
theta_s3 = params_s3['theta_NY']

print(f"\n  V_eff_EC (S3) = {veff_ec_s3}")
print(f"  V_eff_LC (S3) = {veff_lc_s3}")

d3_s3_theta_eta_V_ec = factor(diff(veff_ec_s3, theta_s3, eta_s3, V_s3))
d3_s3_theta_eta_V_lc = factor(diff(veff_lc_s3, theta_s3, eta_s3, V_s3))
print(f"\n  [S3 theta-cubic] d3V_EC/d(theta)d(eta)d(V) = {d3_s3_theta_eta_V_ec}")
print(f"                   d3V_LC/d(theta)d(eta)d(V) = {d3_s3_theta_eta_V_lc}")
diff_theta_s3 = factor(d3_s3_theta_eta_V_ec - d3_s3_theta_eta_V_lc)
print(f"                   EC - LC = {diff_theta_s3}")
pass_theta_s3 = diff_theta_s3.is_zero is True
print(f"  S3 theta-cubic preserved in EC: {'PASS' if pass_theta_s3 else 'FAIL'}")

d3_s3_alpha_V_eta = factor(diff(veff_ec_s3, alpha_s3, V_s3, eta_s3))
alpha_cubic_nonzero_s3 = d3_s3_alpha_V_eta.is_zero is not True
print(f"\n  [S3 alpha-cubic] d3V_EC/d(alpha)d(V)d(eta) = {d3_s3_alpha_V_eta}")
print(f"  S3 alpha-cubic is nonzero: {'PASS' if alpha_cubic_nonzero_s3 else 'FAIL'}")

d3_s3_r_eta_V_ec = factor(diff(veff_ec_s3, r_s3, eta_s3, V_s3))
d3_s3_r_eta_V_lc = factor(diff(veff_lc_s3, r_s3, eta_s3, V_s3))
print(f"\n  [S3 r-eta-V cubic] d3V_EC/d(r)d(eta)d(V) = {d3_s3_r_eta_V_ec}")
print(f"                     d3V_LC/d(r)d(eta)d(V) = {d3_s3_r_eta_V_lc}")
print(f"                     EC - LC = {factor(d3_s3_r_eta_V_ec - d3_s3_r_eta_V_lc)}")

# ── Part 2: T3 cubic vertex analysis ─────────────────────────────────────────
print()
print("-" * 70)
print("Part 2: T3 cubic vertex analysis")
print("-" * 70)
print()

veff_ec_t3, veff_lc_t3, params_t3 = get_veff(TopologyType.T3, "T3")

r_t3   = params_t3['r']
eta_t3 = params_t3['eta']
V_t3   = params_t3['V']
alpha  = params_t3['alpha']
theta  = params_t3['theta_NY']

print(f"\n  V_eff_EC (T3) = {veff_ec_t3}")
print(f"  V_eff_LC (T3) = {veff_lc_t3}")

# theta-cubic: d3V/d(theta)d(eta)d(V)
d3_t3_theta_eta_V_ec = factor(diff(veff_ec_t3, theta, eta_t3, V_t3))
d3_t3_theta_eta_V_lc = factor(diff(veff_lc_t3, theta, eta_t3, V_t3))
print(f"\n  [theta-cubic] d3V_EC/d(theta)d(eta)d(V) = {d3_t3_theta_eta_V_ec}")
print(f"                d3V_LC/d(theta)d(eta)d(V) = {d3_t3_theta_eta_V_lc}")
diff_theta = factor(d3_t3_theta_eta_V_ec - d3_t3_theta_eta_V_lc)
print(f"                EC - LC = {diff_theta}")
pass_theta_preserved = diff_theta.is_zero is True
print(f"  theta-cubic preserved in EC: {'PASS' if pass_theta_preserved else 'FAIL'}")

# alpha-cubic: d3V/d(alpha)d(V)d(eta)
d3_t3_alpha_V_eta = factor(diff(veff_ec_t3, alpha, V_t3, eta_t3))
print(f"\n  [alpha-cubic] d3V_EC/d(alpha)d(V)d(eta) = {d3_t3_alpha_V_eta}")
alpha_cubic_nonzero = d3_t3_alpha_V_eta.is_zero is not True
print(f"  alpha-cubic is nonzero: {'PASS (new EC cubic coefficient)' if alpha_cubic_nonzero else 'FAIL'}")

# r-eta-V mixed cubic
d3_t3_r_eta_V_ec = factor(diff(veff_ec_t3, r_t3, eta_t3, V_t3))
d3_t3_r_eta_V_lc = factor(diff(veff_lc_t3, r_t3, eta_t3, V_t3))
print(f"\n  [r-eta-V cubic] d3V_EC/d(r)d(eta)d(V) = {d3_t3_r_eta_V_ec}")
print(f"                  d3V_LC/d(r)d(eta)d(V) = {d3_t3_r_eta_V_lc}")
print(f"                  EC - LC = {factor(d3_t3_r_eta_V_ec - d3_t3_r_eta_V_lc)}")

# ── Part 3: Pontryagin density alpha-independence ────────────────────────────
print()
print("-" * 70)
print("Part 3: Pontryagin density P alpha-independence")
print("-" * 70)
print()
print("  P(r, eta, V) is computed from MIXING fiber mode cubic coefficient.")
print("  Expected: dP/dalpha = 0 (EC-Weyl does not modify Pontryagin structure)")
print()

cfg_mix = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.MX,
    fiber_mode=FiberMode.MIXING,
    ny_variant=NyVariant.FULL,
)
eng_mix = UnifiedEngine(cfg_mix)
eng_mix.run()

params_mix = eng_mix.data['params']
delta0 = params_mix.get('delta0')
delta1 = params_mix.get('delta1')
delta2 = params_mix.get('delta2')
alpha_mix = params_mix.get('alpha')
V_eff_mix = eng_mix.data['potential']

if delta0 is not None and delta1 is not None and delta2 is not None:
    P_coeff = cancel(
        diff(diff(diff(V_eff_mix, delta0), delta1), delta2).subs(
            [(delta0, 0), (delta1, 0), (delta2, 0)]
        )
    )
    print(f"  P_coeff (cubic coefficient) = {P_coeff}")

    if alpha_mix is not None:
        dP_dalpha = cancel(diff(P_coeff, alpha_mix))
        pass_P_alpha = dP_dalpha.is_zero is True
        print(f"  dP/dalpha = {dP_dalpha}")
        print(f"  P is alpha-independent: {'PASS' if pass_P_alpha else 'FAIL'}")
    else:
        print("  (alpha not in MIXING mode params -- skipping alpha-independence check)")
        pass_P_alpha = None
else:
    print("  MIXING mode delta symbols not available")
    pass_P_alpha = None

# ── Part 4: Nil3 cubic vertex comparison ─────────────────────────────────────
print()
print("-" * 70)
print("Part 4: Nil3 cubic vertex comparison")
print("-" * 70)
print()

veff_ec_nil, veff_lc_nil, params_nil = get_veff(TopologyType.NIL3, "Nil3")

r_nil   = params_nil['r']
eta_nil = params_nil['eta']
V_nil   = params_nil['V']
alpha_nil = params_nil['alpha']
theta_nil = params_nil['theta_NY']

# theta-cubic for Nil3
d3_nil_theta_ec = factor(diff(veff_ec_nil, theta_nil, eta_nil, V_nil))
d3_nil_theta_lc = factor(diff(veff_lc_nil, theta_nil, eta_nil, V_nil))
print(f"  [Nil3 theta-cubic] d3V_EC/d(theta)d(eta)d(V) = {d3_nil_theta_ec}")
print(f"                     d3V_LC/d(theta)d(eta)d(V) = {d3_nil_theta_lc}")
print(f"                     EC - LC = {factor(d3_nil_theta_ec - d3_nil_theta_lc)}")

# alpha-cubic for Nil3
d3_nil_alpha = factor(diff(veff_ec_nil, alpha_nil, V_nil, eta_nil))
print(f"\n  [Nil3 alpha-cubic] d3V_EC/d(alpha)d(V)d(eta) = {d3_nil_alpha}")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("EC Cubic Vertex Summary")
print("=" * 70)
print()
print("  [S3 cubic channels]")
print(f"    theta-cubic (NY channel) preserved: {'PASS' if pass_theta_s3 else 'FAIL'}")
print(f"    alpha-cubic (EC-specific) != 0:      {'PASS' if alpha_cubic_nonzero_s3 else 'FAIL'}")
print()
print("  [T3 cubic vertices]")
print(f"    theta-cubic (NY channel) preserved: {'PASS' if pass_theta_preserved else 'FAIL'}")
print(f"    alpha-cubic (EC-specific) != 0:     {'PASS' if alpha_cubic_nonzero else 'FAIL'}")
print(f"    P alpha-independence:               {'PASS' if pass_P_alpha is True else '? (check)' if pass_P_alpha is None else 'FAIL'}")
print()
print("  [Nil3 cubic channels]")
print(f"    theta-cubic (NY channel) preserved: {'PASS' if factor(d3_nil_theta_ec - d3_nil_theta_lc).is_zero is True else 'FAIL'}")
print(f"    alpha-cubic (EC-specific) != 0:     {'PASS' if d3_nil_alpha.is_zero is not True else 'FAIL'}")
print()
print("  Structural results from computation:")
if alpha_cubic_nonzero_s3:
    print("    1. In S3, EC-Weyl adds a homogeneous EFT cubic coefficient: d3V/d(alpha)d(V)d(eta) != 0")
else:
    print("    1. S3 alpha-cubic = 0  (unexpected: check EC coupling)")
if pass_theta_preserved:
    print("    2. The NY cubic channel is unchanged by EC in T3, and likewise in S3/Nil3")
else:
    print("    2. theta-cubic differs between EC and LC  (unexpected)")
if pass_P_alpha is True:
    print("    3. The MIXING-mode cubic coefficient used for P is alpha-independent")
elif pass_P_alpha is None:
    print("    3. P alpha-independence: not verified (check MIXING mode)")
else:
    print("    3. dP/dalpha != 0  (unexpected: EC modifies Pontryagin)")

teardown_log()
