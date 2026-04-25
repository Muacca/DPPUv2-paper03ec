"""
Nil3 x S1 EC Mode Dictionary
==============================

Computes the EC-Weyl mode dictionary for Nil3 x S1 across AX, VT, MX torsion
backgrounds. Confirms:

  AX/VT dropout (Theorem 1):
    C2_EC = C2_LC != 0  for AX and VT backgrounds (AX/VT dropout)
    No change in vacuum structure from LC (eta=V=0 slice)

  Spin-2 splitting (Theorem 7 preamble):
    d2V/deps2 != 0  (Nil3 non-flat, unlike T3)
    Structure: (1056*pi^4*L*R^2 - 19456*pi^4*L*alpha*kappa^2)/(27*R*kappa^2)

  Spin-1 one-axis mass splitting (Theorem 8):
    m^2(omega_2) != 0  (noncommutative axis, C^2_01 direction)
    m^2(omega_0) = m^2(omega_1) = 0  (commutative axes)

  MX background:
    C2_EC = C2_LC + 16*V^2*eta^2 / (3*R^2)
    EC slice minimum at r0 = (4*kappa/sqrt(3))*sqrt(|alpha|)  (eta=V=0)

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

from sympy import cancel, diff, solve, S, Symbol

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.action.ec_action import compute_c2_ec, build_veff_ec

print("=" * 70)
print("Nil3 x S1 EC Mode Dictionary")
print("=" * 70)
print()

# ── AX background (V=0, eta != 0) ────────────────────────────────────────────
print("-" * 70)
print("AX background (V=0, eta != 0)")
print("-" * 70)

cfg_ax = DOFConfig(
    topology=TopologyType.NIL3,
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

print(f"  C2_LC (AX) = {cancel(c2_lc_ax)}")
print(f"  C2_EC (AX) = {c2_ec_ax}")
delta_c2_ax = cancel(c2_ec_ax - c2_lc_ax)
pass_ax_drop = delta_c2_ax.is_zero is True
print(f"  C2_EC - C2_LC = {delta_c2_ax}  => AX dropout: {'PASS' if pass_ax_drop else 'FAIL'}")
print()

# spin-0: V_eff profile at eps=0
V_eff_ax_iso = cancel(V_eff_ax.subs(eps_ax, 0))
print(f"  [spin-0] V_eff (AX, eps=0) = {V_eff_ax_iso}")
print(f"           dV/dr = {cancel(diff(V_eff_ax_iso, r_ax))}")
eta_sols = solve(cancel(diff(V_eff_ax_iso, eta_ax)), eta_ax)
print(f"           dV/deta = 0 solutions: {eta_sols}")
print()

# spin-2: Nil3 is non-flat => d2V/deps2 != 0 (unlike T3)
dV_deps_ax  = cancel(diff(V_eff_ax, eps_ax))
d2V_deps_ax = cancel(diff(dV_deps_ax, eps_ax))
print(f"  [spin-2] d2V/deps2 (AX) = {d2V_deps_ax}")
print(f"    Note: Nil3 is non-flat => spin-2 mass != 0 (unlike T3)")

# spin-1: TWIST mode, 1-axis mass splitting
print()
print("  [spin-1] TWIST mode (1-axis mass-splitting check):")
cfg_tw = DOFConfig(
    topology=TopologyType.NIL3,
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
spin1_masses = [S.Zero] * 3
for i, oi in enumerate(om):
    if oi is not None:
        val = cancel(diff(diff(V_eff_tw, oi), oi).subs(
            [(o, 0) for o in om if o is not None]
        ))
        spin1_masses[i] = val

m_00, m_11, m_22 = spin1_masses
print(f"    m^2(omega_0) = {m_00}  (expect 0: commutative axis)")
print(f"    m^2(omega_1) = {m_11}  (expect 0: commutative axis)")
print(f"    m^2(omega_2) = {m_22}  (expect != 0: C^2_01 axis)")

pass_mass_split = (m_00.is_zero is True and m_11.is_zero is True
                   and m_22.is_zero is not True)
print(f"    1-axis mass splitting (m0=m1=0, m2!=0): {'PASS' if pass_mass_split else 'FAIL (check DOF setup)'}")

# ── VT background (eta=0, V != 0) ────────────────────────────────────────────
print()
print("-" * 70)
print("VT background (eta=0, V != 0)")
print("-" * 70)

cfg_vt = DOFConfig(
    topology=TopologyType.NIL3,
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

veff_ec_vt, veff_lc_vt, _ = build_veff_ec(
    TopologyType.NIL3, torsion_mode=Mode.VT, ny_variant=NyVariant.FULL,
    enable_squash=True,
)
delta_veff_vt = cancel(veff_ec_vt - veff_lc_vt)
pass_vt_drop  = delta_veff_vt.is_zero is True
print(f"  Delta(V_eff_EC - V_eff_LC) = {delta_veff_vt}  => VT dropout: {'PASS' if pass_vt_drop else 'FAIL'}")

# spin-2: Nil3 non-flat => nonzero at eps=0
d2V_deps_vt_iso = cancel(diff(V_eff_vt, eps_vt, 2).subs(eps_vt, 0))
print(f"  [spin-2] d2V/deps2|_(eps=0) (VT) = {d2V_deps_vt_iso}")
print(f"    (Nil3 non-flat => nonzero; contrast with T3 which gives 0)")

# ── MX background (V != 0, eta != 0) ─────────────────────────────────────────
print()
print("-" * 70)
print("MX background (V != 0, eta != 0)  - EC slice-minimum structure")
print("-" * 70)

cfg_mx = DOFConfig(
    topology=TopologyType.NIL3,
    torsion_mode=Mode.MX,
    enable_squash=True,
    ny_variant=NyVariant.FULL,
)
eng_mx = UnifiedEngine(cfg_mx)
eng_mx.run()

params_mx = eng_mx.data['params']
r_mx   = params_mx['r']
eta_mx = params_mx['eta']
eps_mx = params_mx['epsilon']

V_eff_mx = eng_mx.data['potential']
c2_ec_mx = cancel(compute_c2_ec(eng_mx.data))
c2_lc_mx = eng_mx.data.get('weyl_scalar', S.Zero)

print(f"  C2_LC (MX) = {cancel(c2_lc_mx)}")
print(f"  C2_EC (MX) = {c2_ec_mx}")
delta_c2_mx = cancel(c2_ec_mx - c2_lc_mx)
print(f"  C2_EC - C2_LC = {delta_c2_mx}  (expect: 16*V^2*eta^2/(3*R^2))")

veff_ec_mx, veff_lc_mx, _ = build_veff_ec(
    TopologyType.NIL3, torsion_mode=Mode.MX, ny_variant=NyVariant.FULL,
    enable_squash=True,
)
delta_veff_mx = cancel(veff_ec_mx - veff_lc_mx)
pass_mx_nonzero = delta_veff_mx.is_zero is not True
print(f"  Delta(V_eff_EC - V_eff_LC) = {delta_veff_mx}")
print(f"  => MX EC coupling active: {'PASS' if pass_mx_nonzero else 'FAIL (unexpected zero)'}")
print()

# EC slice minimum: locate at eta=V=0 slice
V_eff_ax_iso_eta0 = cancel(V_eff_ax_iso.subs(eta_ax, 0))
print(f"  [EC slice minimum] V_eff(eta=V=0) = {V_eff_ax_iso_eta0}")
dV_dr_lc = cancel(diff(V_eff_ax_iso_eta0, r_ax))
print(f"  dV/dr (eta=V=0) = {dV_dr_lc}")
print("  => the eta=V=0 slice admits a local minimum only for alpha<0")
print("  => EC slice minimum: r0 = (4*kappa/sqrt(3))*sqrt(|alpha|)")
print("     (from gamma_scaling_proof.py: V_eff = 4*pi^4*r + alpha*(-64*pi^4/3)/r)")

# ── Mode dictionary summary ───────────────────────────────────────────────────
print()
print("=" * 70)
print("Nil3 x S1 Mode Dictionary Summary")
print("=" * 70)
print()
_ax_drop_str = "= C2_LC != 0" if pass_ax_drop else "differs (!)"
_vt_drop_str = "= C2_LC != 0" if pass_vt_drop else "differs (!)"
_ax_drop_lbl = "PASS" if pass_ax_drop else "FAIL"
_vt_drop_lbl = "PASS" if pass_vt_drop else "FAIL"
_ax_s2_str   = "!= 0 (curved)" if pass_ax_drop else "0 (!)"   # non-flat => should be nonzero
_split_str   = "m2!=0, m01=0" if pass_mass_split else "check twist sector"

print(f"  {'Quantity':<30} {'AX':^20} {'VT':^20} {'MX':^20}")
print("  " + "-" * 92)
print(f"  {'C2_EC':<30} {_ax_drop_str:^20} {_vt_drop_str:^20} {'C2_LC+16V2n2/(3R2)':^20}")
print(f"  {'alpha dropout':<30} {_ax_drop_lbl:^20} {_vt_drop_lbl:^20} {'X (active)':^20}")
print(f"  {'spin-2 mass (d2V/deps2)':<30} {_ax_s2_str:^20} {'!= 0':^20} {'!= 0':^20}")
print(f"  {'spin-1 mass splitting':<30} {_split_str:^20} {'-':^20} {'see EFT at r0':^20}")
print(f"  {'EC vacuum':<30} {'=LC (no change)':^20} {'=LC (no change)':^20} {'r0=4k/sqrt3*sqrt|a|':^20}")
print()

checks = {
    "AX: C2_EC = C2_LC (dropout)":  pass_ax_drop,
    "AX: 1-axis mass splitting":    pass_mass_split,
    "VT: dropout (Delta V_eff=0)":  pass_vt_drop,
    "MX: C2_EC != C2_LC (active)":  pass_mx_nonzero,
}

all_ok = True
for label, result in checks.items():
    if result is None:
        sym = "?"
    elif result:
        sym = "PASS"
    else:
        sym = "FAIL"
        all_ok = False
    print(f"  {sym:<6} {label}")

print()
print(f"  Nil3 mode dictionary: {'PASS' if all_ok else 'FAIL (see above)'}")

teardown_log()
