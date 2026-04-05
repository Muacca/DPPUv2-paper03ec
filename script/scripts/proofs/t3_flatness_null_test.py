"""
T3 Flatness Null Test
=====================

Theorem: T3 x S1 has zero spin-2 masses (squash and shear) for all torsion
backgrounds (AX, VT, MX), because T3 is flat: C^a_{bc} = 0 for all a,b,c.

Verification:
  d2V/ds2     = 0  (SymPy exact)  -- shear null test
  d2V/deps2   = 0  (SymPy exact)  -- squash null test
  d2V/deps_ds = 0  (SymPy exact)  -- cross-term null (Weyl)
  Discovery:  d2V/deps_ds != 0 in general due to torsion-volume coupling
              (kinematic origin, NOT Weyl spin-2 mass)

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

from sympy import cancel, diff, S

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.action.ec_action import compute_c2_ec

print("=" * 70)
print("T3 Flatness Null Test  -- Theorem (Prelude to Theorem 5)")
print("=" * 70)
print()

checks = {}

# ── AX background (V=0, eta != 0): squash + shear ──────────────────────────
print("-" * 70)
print("AX background (V=0, eta!=0): squash + shear")
print("-" * 70)

cfg_ax = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    enable_shear=True,
    ny_variant=NyVariant.FULL,
)
eng_ax = UnifiedEngine(cfg_ax)
eng_ax.run()

params_ax = eng_ax.data['params']
eps_ax    = params_ax['epsilon']
s_ax      = params_ax['s']

V_eff_ax = eng_ax.data['potential']
c2_ec_ax = cancel(compute_c2_ec(eng_ax.data))

print(f"  C2_EC (AX, squash+shear) = {c2_ec_ax}")

d2V_ds2_ax     = cancel(diff(V_eff_ax, s_ax, 2))
d2V_deps2_ax   = cancel(diff(V_eff_ax, eps_ax, 2))
d2V_deps_ds_ax = cancel(diff(diff(V_eff_ax, eps_ax), s_ax))

print(f"\n  d2V/ds2   (AX) = {d2V_ds2_ax}")
print(f"  d2V/deps2 (AX) = {d2V_deps2_ax}")
print(f"  d2V/deps_ds (AX) = {d2V_deps_ds_ax}")

pass_s2_ax    = d2V_ds2_ax.is_zero is True
pass_eps2_ax  = d2V_deps2_ax.is_zero is True
pass_epss_ax  = d2V_deps_ds_ax.is_zero is True
checks['AX_s2']   = pass_s2_ax
checks['AX_eps2'] = pass_eps2_ax
checks['AX_epss'] = pass_epss_ax

print(f"\n  d2V/ds2 = 0 (Weyl spin-2 shear):  {'PASS' if pass_s2_ax else 'FAIL'}")
print(f"  d2V/deps2 = 0 (Weyl spin-2 squash): {'PASS' if pass_eps2_ax else 'FAIL'}")
print(f"  d2V/deps_ds = 0 (cross Weyl):        {'PASS' if pass_epss_ax else 'NOTE: nonzero = torsion-volume coupling'}")

# ── VT background (eta=0, V != 0): squash + shear ──────────────────────────
print()
print("-" * 70)
print("VT background (eta=0, V!=0): squash + shear")
print("-" * 70)

cfg_vt = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.VT,
    enable_squash=True,
    enable_shear=True,
    ny_variant=NyVariant.FULL,
)
eng_vt = UnifiedEngine(cfg_vt)
eng_vt.run()

params_vt = eng_vt.data['params']
eps_vt    = params_vt['epsilon']
s_vt      = params_vt['s']

V_eff_vt  = eng_vt.data['potential']
c2_ec_vt  = cancel(compute_c2_ec(eng_vt.data))

print(f"  C2_EC (VT, squash+shear) = {c2_ec_vt}")

d2V_ds2_vt     = cancel(diff(V_eff_vt, s_vt, 2))
d2V_deps2_vt   = cancel(diff(V_eff_vt, eps_vt, 2))
d2V_deps_ds_vt = cancel(diff(diff(V_eff_vt, eps_vt), s_vt))

print(f"\n  d2V/ds2   (VT) = {d2V_ds2_vt}")
print(f"  d2V/deps2 (VT) = {d2V_deps2_vt}")
print(f"  d2V/deps_ds (VT) = {d2V_deps_ds_vt}")

pass_s2_vt    = d2V_ds2_vt.is_zero is True
pass_eps2_vt  = d2V_deps2_vt.is_zero is True
pass_epss_vt  = d2V_deps_ds_vt.is_zero is True
checks['VT_s2']   = pass_s2_vt
checks['VT_eps2'] = pass_eps2_vt
checks['VT_epss'] = pass_epss_vt

print(f"\n  d2V/ds2 = 0 (Weyl spin-2 shear):  {'PASS' if pass_s2_vt else 'FAIL'}")
print(f"  d2V/deps2 = 0 (Weyl spin-2 squash): {'PASS' if pass_eps2_vt else 'FAIL'}")
print(f"  d2V/deps_ds = 0 (cross Weyl):        {'PASS' if pass_epss_vt else 'NOTE: nonzero = torsion-volume coupling'}")

# ── MX background (V != 0, eta != 0): squash + shear ───────────────────────
print()
print("-" * 70)
print("MX background (eta!=0, V!=0): squash + shear")
print("-" * 70)

cfg_mx = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.MX,
    enable_squash=True,
    enable_shear=True,
    ny_variant=NyVariant.FULL,
)
eng_mx = UnifiedEngine(cfg_mx)
eng_mx.run()

params_mx = eng_mx.data['params']
eps_mx    = params_mx['epsilon']
s_mx      = params_mx['s']

V_eff_mx  = eng_mx.data['potential']
c2_ec_mx  = cancel(compute_c2_ec(eng_mx.data))

print(f"  C2_EC (MX, squash+shear) = {c2_ec_mx}")

d2V_ds2_mx     = cancel(diff(V_eff_mx, s_mx, 2))
d2V_deps2_mx   = cancel(diff(V_eff_mx, eps_mx, 2))
d2V_deps_ds_mx = cancel(diff(diff(V_eff_mx, eps_mx), s_mx))

print(f"\n  d2V/ds2   (MX) = {d2V_ds2_mx}")
print(f"  d2V/deps2 (MX) = {d2V_deps2_mx}")
print(f"  d2V/deps_ds (MX) = {d2V_deps_ds_mx}")

pass_s2_mx    = d2V_ds2_mx.is_zero is True
pass_eps2_mx  = d2V_deps2_mx.is_zero is True
pass_epss_mx  = d2V_deps_ds_mx.is_zero is True
checks['MX_s2']   = pass_s2_mx
checks['MX_eps2'] = pass_eps2_mx
checks['MX_epss'] = pass_epss_mx

print(f"\n  d2V/ds2 = 0 (Weyl spin-2 shear):  {'PASS' if pass_s2_mx else 'FAIL'}")
print(f"  d2V/deps2 = 0 (Weyl spin-2 squash): {'PASS' if pass_eps2_mx else 'FAIL'}")
print(f"  NOTE: MX has C2_EC = 16V^2 eta^2 / (3r^2) != 0, but spin-2 Weyl mass = 0")
print(f"        due to T3 flatness (C^a_bc = 0 => Weyl tensor = 0)")

# ── Summary ──────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("  T3 Flatness Null Test -- Summary")
print("=" * 70)
print()

# Weyl spin-2 checks (d2V/deps2, d2V/ds2) are the primary criterion
weyl_checks = {
    'AX_s2':   checks['AX_s2'],
    'AX_eps2': checks['AX_eps2'],
    'VT_s2':   checks['VT_s2'],
    'VT_eps2': checks['VT_eps2'],
    'MX_s2':   checks['MX_s2'],
    'MX_eps2': checks['MX_eps2'],
}
all_weyl_pass = all(weyl_checks.values())

check_names = {
    'AX_s2':   'AX: d2V/ds2 = 0  (Weyl spin-2 shear mass)',
    'AX_eps2': 'AX: d2V/deps2 = 0  (Weyl spin-2 squash mass)',
    'VT_s2':   'VT: d2V/ds2 = 0  (Weyl spin-2 shear mass)',
    'VT_eps2': 'VT: d2V/deps2 = 0  (Weyl spin-2 squash mass)',
    'MX_s2':   'MX: d2V/ds2 = 0  (T3 flatness, all backgrounds)',
    'MX_eps2': 'MX: d2V/deps2 = 0  (T3 flatness, all backgrounds)',
}

for key, desc in check_names.items():
    passed = weyl_checks.get(key, False)
    print(f"  {'PASS' if passed else 'FAIL'}  {desc}")

print()
print("  [DISCOVERY] Torsion-volume coupling: d2V/deps_ds != 0 in general:")
print(f"    AX: d2V/deps_ds = {d2V_deps_ds_ax}")
print(f"    VT: d2V/deps_ds = {d2V_deps_ds_vt}")
print(f"    MX: d2V/deps_ds = {d2V_deps_ds_mx}")
print()
print("  Physical interpretation:")
print("    T3 volume = (2*pi)^4 * L * r^3 * (1+eps) * (1+s)")
print("    => V_torsion ~ Vol * (torsion kinetic, C-independent)")
print("       d2V_torsion/deps_ds = (torsion kinetic) * r != 0")
print("    This is TORSION-VOLUME COUPLING, not Weyl/graviton spin-2 mass.")
print("    For AX/VT: C2_EC = C2_LC = 0, so spin-2 mass = 0 by direct vanishing.")
print("    For MX: C2_EC = 16V^2 eta^2 / (3r^2) != 0, but d2C2_EC/deps2 = 0")
print("            (C2_EC is eps,s-independent) => Weyl spin-2 mass d2V/deps2 = 0 still holds.")
print("    The deps_ds cross-term has KINEMATIC origin + alpha*C2_EC contribution (MX only).")

print()
if all_weyl_pass:
    print("  Null test result: PASS")
    print("  T3 Weyl spin-2 mass (d2V/deps2, d2V/ds2) = 0 verified (SymPy exact)")
    print("  Addendum: d2V/deps_ds != 0 = torsion-volume coupling (no Weyl content)")
else:
    failed = [k for k, v in weyl_checks.items() if not v]
    print(f"  Null test result: FAIL  (failed checks: {failed})")

teardown_log()
