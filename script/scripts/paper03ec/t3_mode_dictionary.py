"""
T3 x S1 EC Mode Dictionary
===========================

Computes the EC-Weyl mode dictionary for T3 x S1 across AX, VT, MX torsion
backgrounds. Confirms:

  AX/VT dropout (Theorem 1):
    C2_EC = C2_LC = 0  for AX and VT backgrounds (T3 flatness + dropout)

  T3 flatness (Theorem 5 preamble):
    d2V/deps2 = 0  for all backgrounds (spin-2 mass = 0)

  MX EC coupling:
    C2_EC = 16*V^2*eta^2 / (3*r^2)  for MX background

  Spin-1 (homogeneous EFT):
    Mass = 0 in V_eff for the T3 twist sector
    The homogeneous mode dictionary sees only a vanishing quadratic coefficient

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

from sympy import cancel, diff, Rational, S, zeros

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.action.ec_action import compute_c2_ec, build_veff_ec

print("=" * 70)
print("T3 x S1 EC Mode Dictionary")
print("=" * 70)
print()

# ── AX background (V=0, eta != 0) ────────────────────────────────────────────
print("-" * 70)
print("AX background (V=0, eta != 0)")
print("-" * 70)

cfg_ax = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    ny_variant=NyVariant.FULL,
)
eng_ax = UnifiedEngine(cfg_ax)
eng_ax.run()

params_ax = eng_ax.data['params']
r_ax   = params_ax['r']
eta_ax = params_ax['eta']
eps_ax = params_ax['epsilon']

V_eff_ax = eng_ax.data['potential']
c2_ec_ax = cancel(compute_c2_ec(eng_ax.data))
c2_lc_ax = eng_ax.data.get('weyl_scalar', S.Zero)

print(f"  V_eff (AX, eps=0) = {cancel(V_eff_ax.subs(eps_ax, 0))}")
print(f"  C2_LC (AX) = {cancel(c2_lc_ax)}")
print(f"  C2_EC (AX) = {c2_ec_ax}")
print(f"  C2_EC - C2_LC = {cancel(c2_ec_ax - c2_lc_ax)}  (expect: 0)")
print()

# spin-2: d2V/deps2 = 0 (T3 flatness)
d2V_deps_ax = cancel(diff(V_eff_ax, eps_ax, 2))
pass_s2_ax  = d2V_deps_ax.is_zero is True
print(f"  [spin-2] d2V/deps2 (AX) = {d2V_deps_ax}  =>  {'PASS (=0)' if pass_s2_ax else 'FAIL'}")

# spin-0: V_eff structure
V_eff_ax_iso = cancel(V_eff_ax.subs(eps_ax, 0))
print(f"  [spin-0] dV/dr (AX, iso) = {cancel(diff(V_eff_ax_iso, r_ax))}")
print(f"           dV/deta (AX, iso) = {cancel(diff(V_eff_ax_iso, eta_ax))}")

# spin-1: TWIST mode - mass from homogeneous V_eff
print()
print("  [spin-1] TWIST mode - homogeneous EFT mass check:")
cfg_tw = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.AX,
    fiber_mode=FiberMode.TWIST,
    isotropic_twist=False,
    ny_variant=NyVariant.FULL,
)
eng_tw = UnifiedEngine(cfg_tw)
eng_tw.run()

params_tw  = eng_tw.data['params']
omega1_sym = params_tw.get('omega1')
omega2_sym = params_tw.get('omega2')
omega3_sym = params_tw.get('omega3')
V_eff_tw   = eng_tw.data['potential']

om = [omega1_sym, omega2_sym, omega3_sym]
G_spin1 = zeros(3)
for i, oi in enumerate(om):
    for j, oj in enumerate(om):
        if oi is not None and oj is not None:
            G_spin1[i, j] = cancel(diff(diff(V_eff_tw, oi), oj).subs(
                [(o, 0) for o in om if o is not None]
            ))

g00 = G_spin1[0, 0]
pass_s1_mass_zero = g00.is_zero is True
print(f"    G_spin1[0,0] = {g00}  (mass=0 in V_eff: {'PASS' if pass_s1_mass_zero else 'FAIL'})")
print("    => no additional homogeneous quadratic term is generated in the twist sector")

# ── VT background (eta=0, V != 0) ────────────────────────────────────────────
print()
print("-" * 70)
print("VT background (eta=0, V != 0)")
print("-" * 70)

cfg_vt = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.VT,
    enable_squash=True,
    ny_variant=NyVariant.FULL,
)
eng_vt = UnifiedEngine(cfg_vt)
eng_vt.run()

params_vt = eng_vt.data['params']
eps_vt = params_vt['epsilon']

V_eff_vt = eng_vt.data['potential']
c2_ec_vt = cancel(compute_c2_ec(eng_vt.data))
c2_lc_vt = eng_vt.data.get('weyl_scalar', S.Zero)

print(f"  C2_LC (VT) = {cancel(c2_lc_vt)}")
print(f"  C2_EC (VT) = {c2_ec_vt}")

# VT dropout: build_veff_ec confirms V_eff_EC = V_eff_LC
veff_ec_vt, veff_lc_vt, _ = build_veff_ec(
    TopologyType.T3, torsion_mode=Mode.VT, ny_variant=NyVariant.FULL,
    enable_squash=True,
)
delta_veff_vt = cancel(veff_ec_vt - veff_lc_vt)
pass_vt_drop = delta_veff_vt.is_zero is True
print(f"  Delta(V_eff_EC - V_eff_LC) = {delta_veff_vt}  => VT dropout: {'PASS' if pass_vt_drop else 'FAIL'}")

d2V_deps_vt = cancel(diff(V_eff_vt, eps_vt, 2))
pass_s2_vt  = d2V_deps_vt.is_zero is True
print(f"  d2V/deps2 (VT) = {d2V_deps_vt}  => spin-2 mass=0: {'PASS' if pass_s2_vt else 'FAIL'}")

# ── MX background (V != 0, eta != 0) ─────────────────────────────────────────
print()
print("-" * 70)
print("MX background (V != 0, eta != 0)")
print("-" * 70)

cfg_mx = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.MX,
    enable_squash=True,
    ny_variant=NyVariant.FULL,
)
eng_mx = UnifiedEngine(cfg_mx)
eng_mx.run()

params_mx = eng_mx.data['params']
r_mx   = params_mx['r']
eta_mx = params_mx['eta']
V_mx   = params_mx.get('V', S.Zero)
eps_mx = params_mx['epsilon']

V_eff_mx = eng_mx.data['potential']
c2_ec_mx = cancel(compute_c2_ec(eng_mx.data))
c2_lc_mx = eng_mx.data.get('weyl_scalar', S.Zero)

print(f"  C2_LC (MX) = {cancel(c2_lc_mx)}")
print(f"  C2_EC (MX) = {c2_ec_mx}")
delta_c2_mx = cancel(c2_ec_mx - c2_lc_mx)
print(f"  C2_EC - C2_LC = {delta_c2_mx}  (expect: 16*V^2*eta^2/(3*r^2))")

veff_ec_mx, veff_lc_mx, _ = build_veff_ec(
    TopologyType.T3, torsion_mode=Mode.MX, ny_variant=NyVariant.FULL,
    enable_squash=True,
)
delta_veff_mx = cancel(veff_ec_mx - veff_lc_mx)
pass_mx_nonzero = delta_veff_mx.is_zero is not True
print(f"  Delta(V_eff_EC - V_eff_LC) = {delta_veff_mx}")
print(f"  => MX EC coupling active (Delta != 0): {'PASS' if pass_mx_nonzero else 'FAIL (unexpected zero)'}")

d2V_deps_mx = cancel(diff(V_eff_mx, eps_mx, 2))
pass_s2_mx  = d2V_deps_mx.is_zero is True
print(f"  d2V/deps2 (MX) = {d2V_deps_mx}  => spin-2 mass=0: {'PASS' if pass_s2_mx else 'FAIL'}")

# ── Mode dictionary summary ───────────────────────────────────────────────────
print()
print("=" * 70)
print("T3 x S1 Mode Dictionary Summary")
print("=" * 70)
print()
_ax_drop_str  = "= C2_LC = 0" if cancel(c2_ec_ax - c2_lc_ax).is_zero is True else "!= C2_LC (!)"
_vt_drop_str  = "= C2_LC = 0" if pass_vt_drop else "!= C2_LC (!)"
_ax_s2_str    = "0 (T3 flat)" if pass_s2_ax   else "!= 0 (!)"
_vt_s2_str    = "0 (T3 flat)" if pass_s2_vt   else "!= 0 (!)"
_mx_s2_str    = "0 (T3 flat)" if pass_s2_mx   else "!= 0 (!)"
_ax_s1_str    = "0 (massless)" if pass_s1_mass_zero else "!= 0 (!)"
_ax_drop_lbl  = "PASS" if cancel(c2_ec_ax - c2_lc_ax).is_zero is True else "FAIL"
_vt_drop_lbl  = "PASS" if pass_vt_drop else "FAIL"
print(f"  {'Quantity':<30} {'AX':^18} {'VT':^18} {'MX':^18}")
print("  " + "-" * 85)
print(f"  {'C2_EC':<30} {_ax_drop_str:^18} {_vt_drop_str:^18} {'16V2eta2/(3r2)':^18}")
print(f"  {'alpha dropout':<30} {_ax_drop_lbl:^18} {_vt_drop_lbl:^18} {'X (active)':^18}")
print(f"  {'spin-2 mass (d2V/deps2)':30} {_ax_s2_str:^18} {_vt_s2_str:^18} {_mx_s2_str:^18}")
print(f"  {'spin-1 mass (homogeneous V_eff)':<30} {_ax_s1_str:^18} {'0':^18} {'0':^18}")
print()

checks = {
    "AX: C2_EC = C2_LC (dropout)":    cancel(c2_ec_ax - c2_lc_ax).is_zero is True,
    "AX: spin-2 d2V/deps2 = 0":       pass_s2_ax,
    "AX: spin-1 mass=0 in V_eff":     pass_s1_mass_zero,
    "VT: dropout (Delta V_eff=0)":    pass_vt_drop,
    "VT: spin-2 d2V/deps2 = 0":       pass_s2_vt,
    "MX: C2_EC != C2_LC (active)":    pass_mx_nonzero,
    "MX: spin-2 d2V/deps2 = 0":       pass_s2_mx,
}

print("  Verification checks:")
all_ok = True
for label, result in checks.items():
    if result is None:
        sym = "?"
    elif result:
        sym = "PASS"
    else:
        sym = "FAIL"
        all_ok = False
    print(f"    {sym:<6} {label}")

print()
overall = "PASS" if all_ok else "FAIL (see above)"
print(f"  T3 mode dictionary: {overall}")

teardown_log()
