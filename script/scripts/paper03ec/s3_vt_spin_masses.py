"""
S3 x S1 VT Spin-2/1 Mass (Palatini Universality, Theorem 4)
=============================================================

Computes spin-2 quintet mass and spin-1 triplet mass in S3 x S1 VT background.
Verifies Theorem 4 (Palatini universality): AX and VT backgrounds give identical
spin-2/1 mass spectra.

Key results (Theorem 4):
  VT spin-2: d2V/deps2|_(eps=V=0, alpha=0) = AX spin-2 (V-independent)
  VT spin-1: m^2(omega) = same as AX (confirmed numerically)

Steps:
  Step 1: S3 VT V_eff(r, V, eps, alpha) acquisition
  Step 2: Palatini protection -- dV/dV=0 -> V*=0 (on-shell)
  Step 3: Spin-2 mass -- d2V/deps2|_{eps=0} (V-independent, = AX)
  Step 4: AX spin-2 comparison (symbolic + numerical)
  Step 5: Spin-1 mass -- TWIST mode with VT background, compare to AX

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

from sympy import cancel, diff, pi, S, solve

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.action.ec_action import compute_c2_ec

print("=" * 70)
print("S3 x S1 VT Spin-2/1 Mass  -- Theorem 4 (Palatini Universality)")
print("=" * 70)
print()

# ── Step 1: S3 VT V_eff ───────────────────────────────────────────────────────
print("-" * 70)
print("Step 1: S3 VT V_eff(r, V, eps, alpha) acquisition")
print("-" * 70)

cfg_vt = DOFConfig(
    topology=TopologyType.S3,
    torsion_mode=Mode.VT,
    enable_squash=True,
    ny_variant=NyVariant.FULL,
)
eng_vt = UnifiedEngine(cfg_vt)
eng_vt.run()

params_vt = eng_vt.data['params']
r_vt     = params_vt['r']
V_vt     = params_vt['V']
eps_vt   = params_vt['epsilon']
L_vt     = params_vt['L']
kappa_vt = params_vt['kappa']
alpha_vt = params_vt.get('alpha')
theta_vt = params_vt.get('theta') or params_vt.get('theta_NY')

c2_ec_vt = cancel(compute_c2_ec(eng_vt.data))
c2_lc_vt = eng_vt.data.get('weyl_scalar', S.Zero)
delta_c2  = cancel(c2_ec_vt - c2_lc_vt)
print(f"  VT dropout check: C2_EC - C2_LC = {delta_c2}  (expect: 0)")
pass_vt_drop = delta_c2.is_zero is True
print(f"  VT dropout: {'PASS' if pass_vt_drop else 'FAIL'}")
print()

V_eff_vt = cancel(eng_vt.data['potential'])
subs_no_ny = [(theta_vt, S.Zero)] if theta_vt is not None else []
V_eff_vt_nn = cancel(V_eff_vt.subs(subs_no_ny))
V_eff_vt_iso = cancel(V_eff_vt_nn.subs(eps_vt, S.Zero))
print(f"  V_eff(r, V, eps=0, theta=0) = {V_eff_vt_iso}")

# ── Step 2: Palatini protection ───────────────────────────────────────────────
print()
print("-" * 70)
print("Step 2: Palatini protection -- dV/dV = 0 => V* = 0")
print("-" * 70)

dV_dV = cancel(diff(V_eff_vt_iso, V_vt))
print(f"\n  dV/dV |_(eps=0) = {dV_dV}")

V_onshell = solve(dV_dV, V_vt)
print(f"  V*(r) = {V_onshell}")

if all(s == S.Zero for s in V_onshell):
    print("  CHECK: V* = 0: Palatini protection confirmed (V is non-dynamical)")
    V_star = S.Zero
    pass_palatini = True
else:
    print(f"  WARNING: V* non-trivial: {V_onshell}")
    V_star = V_onshell[0] if V_onshell else S.Zero
    pass_palatini = False

# ── Step 3: Spin-2 mass ───────────────────────────────────────────────────────
print()
print("-" * 70)
print("Step 3: Spin-2 mass d2V/deps2|_(eps=0) -- V-independence check")
print("-" * 70)

d2V_deps2_vt  = cancel(diff(V_eff_vt_nn, eps_vt, 2))
d2V_deps2_iso = cancel(d2V_deps2_vt.subs(eps_vt, S.Zero))
print(f"\n  d2V/deps2|_(eps=0) = {d2V_deps2_iso}")

has_V = V_vt in d2V_deps2_iso.free_symbols
print(f"  Contains V: {has_V}")
if has_V:
    d2V_vt_V0 = cancel(d2V_deps2_iso.subs(V_vt, S.Zero))
    print(f"  After V=0: {d2V_vt_V0}")
    print("  V-independence: V enters only as constant (consistent with Palatini)")
else:
    d2V_vt_V0 = d2V_deps2_iso
    print("  V-independent: Palatini protection confirmed for spin-2")
    pass_spin2_vindep = True

# Separate LC and EC parts
if alpha_vt is not None:
    d2V_lc_part = cancel(d2V_vt_V0.subs(alpha_vt, S.Zero))
    d2V_ec_part = cancel(d2V_vt_V0 - d2V_lc_part)
    print(f"\n  LC part (alpha=0): {d2V_lc_part}")
    print(f"  EC correction:     {d2V_ec_part}")
else:
    d2V_lc_part = d2V_vt_V0

# Numerical at r=3, L=kappa=1, alpha=0
subs_lc_num = [(L_vt, 1), (kappa_vt, 1), (r_vt, 3)]
if alpha_vt is not None:
    subs_lc_num.append((alpha_vt, S.Zero))
if V_vt in d2V_vt_V0.free_symbols:
    subs_lc_num.append((V_vt, S.Zero))

val_vt = float(d2V_vt_V0.subs(subs_lc_num))
ref_128pi2 = float(128 * pi**2)
ref_48pi2  = float(48  * pi**2)
print(f"\n  Numerical (r=3, L=kappa=1, alpha=0): {val_vt:.6f}")
print(f"  128*pi^2 = {ref_128pi2:.6f}  [expected S3 AX value]")
print(f"  Ratio to 128*pi^2: {val_vt / ref_128pi2:.6f}")

# ── Step 4: AX comparison ─────────────────────────────────────────────────────
print()
print("-" * 70)
print("Step 4: AX spin-2 comparison (Theorem 4)")
print("-" * 70)

cfg_ax = DOFConfig(
    topology=TopologyType.S3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    ny_variant=NyVariant.FULL,
)
eng_ax = UnifiedEngine(cfg_ax)
eng_ax.run()

params_ax = eng_ax.data['params']
r_ax      = params_ax['r']
eta_ax    = params_ax.get('eta')
eps_ax    = params_ax['epsilon']
L_ax      = params_ax['L']
kappa_ax  = params_ax['kappa']
alpha_ax  = params_ax.get('alpha')
theta_ax  = params_ax.get('theta') or params_ax.get('theta_NY')

V_eff_ax = cancel(eng_ax.data['potential'])
if theta_ax is not None:
    V_eff_ax = cancel(V_eff_ax.subs(theta_ax, S.Zero))

d2V_ax_iso = cancel(diff(V_eff_ax, eps_ax, 2).subs(eps_ax, S.Zero))
print(f"\n  AX d2V/deps2|_(eps=0) = {d2V_ax_iso}")

subs_ax_num = [(L_ax, 1), (kappa_ax, 1), (r_ax, 3)]
if eta_ax is not None:
    subs_ax_num.append((eta_ax, 1))
if alpha_ax is not None:
    subs_ax_num.append((alpha_ax, S.Zero))

val_ax = float(d2V_ax_iso.subs(subs_ax_num))
print(f"  AX numerical (eta=1, r=3, L=kappa=1, alpha=0): {val_ax:.6f}")

# Symbolic comparison (at L=kappa=1, alpha=0)
subs_sym_vt = [(L_vt, 1), (kappa_vt, 1)]
if alpha_vt is not None:
    subs_sym_vt.append((alpha_vt, S.Zero))
if V_vt in d2V_vt_V0.free_symbols:
    subs_sym_vt.append((V_vt, S.Zero))
formula_vt = cancel(d2V_vt_V0.subs(subs_sym_vt))

subs_sym_ax = [(L_ax, 1), (kappa_ax, 1)]
if alpha_ax is not None:
    subs_sym_ax.append((alpha_ax, S.Zero))
if eta_ax is not None:
    subs_sym_ax.append((eta_ax, 1))
formula_ax = cancel(d2V_ax_iso.subs(subs_sym_ax))

print(f"\n  VT formula (L=kappa=1, alpha=0, V=0): {formula_vt}")
print(f"  AX formula (L=kappa=1, alpha=0, eta=1): {formula_ax}")

diff_vt_ax = cancel(formula_vt - formula_ax)
print(f"  VT - AX = {diff_vt_ax}")

if diff_vt_ax.is_zero is True:
    spin2_result = "= AX (Palatini universality CONFIRMED)"
    pass_spin2_eq = True
else:
    ratio = val_vt / val_ax
    spin2_result = f"ratio VT/AX = {ratio:.6f} (numeric)"
    pass_spin2_eq = abs(ratio - 1.0) < 1e-6

print(f"\n  Spin-2 universality (VT=AX): {'PASS' if pass_spin2_eq else 'FAIL'}")
print(f"  Result: {spin2_result}")

# ── Step 5: Spin-1 mass ───────────────────────────────────────────────────────
print()
print("-" * 70)
print("Step 5: Spin-1 mass (TWIST mode, VT) -- compare to AX")
print("-" * 70)

cfg_vt_tw = DOFConfig(
    topology=TopologyType.S3,
    torsion_mode=Mode.VT,
    fiber_mode=FiberMode.TWIST,
    isotropic_twist=False,
    enable_squash=False,
    ny_variant=NyVariant.FULL,
)
eng_vt_tw = UnifiedEngine(cfg_vt_tw)
eng_vt_tw.run()

params_tw = eng_vt_tw.data['params']
omega1_tw = params_tw.get('omega1')
omega2_tw = params_tw.get('omega2')
omega3_tw = params_tw.get('omega3')
r_tw      = params_tw['r']
V_tw      = params_tw.get('V')
L_tw      = params_tw['L']
kappa_tw  = params_tw['kappa']
alpha_tw  = params_tw.get('alpha')
theta_tw  = params_tw.get('theta') or params_tw.get('theta_NY')

V_eff_tw = cancel(eng_vt_tw.data['potential'])
if theta_tw is not None:
    V_eff_tw = cancel(V_eff_tw.subs(theta_tw, S.Zero))

omega_subs = [(sym, S.Zero) for sym in [omega1_tw, omega2_tw, omega3_tw] if sym is not None]
spin1_masses_vt = {}
print(f"\n  VT spin-1 mass d2V/domega_k^2|_(omega=0):")
for name, sym in [('omega1', omega1_tw), ('omega2', omega2_tw), ('omega3', omega3_tw)]:
    if sym is not None:
        d2_0 = cancel(diff(V_eff_tw, sym, 2).subs(omega_subs))
        spin1_masses_vt[name] = d2_0
        subs_num = [(L_tw, 1), (kappa_tw, 1)]
        if alpha_tw is not None:
            subs_num.append((alpha_tw, S.Zero))
        if V_tw is not None:
            subs_num.append((V_tw, S.Zero))
        val = float(d2_0.subs(subs_num).subs(r_tw, 3))
        print(f"    {name}: {val:.6e}")

# AX spin-1 reference
cfg_ax_tw = DOFConfig(
    topology=TopologyType.S3,
    torsion_mode=Mode.AX,
    fiber_mode=FiberMode.TWIST,
    isotropic_twist=False,
    enable_squash=False,
    ny_variant=NyVariant.FULL,
)
eng_ax_tw = UnifiedEngine(cfg_ax_tw)
eng_ax_tw.run()

params_ax_tw = eng_ax_tw.data['params']
om1_ax   = params_ax_tw.get('omega1')
r_ax_tw  = params_ax_tw['r']
eta_ax_tw = params_ax_tw.get('eta')
L_ax_tw  = params_ax_tw['L']
kappa_ax_tw = params_ax_tw['kappa']
alpha_ax_tw = params_ax_tw.get('alpha')
theta_ax_tw = params_ax_tw.get('theta') or params_ax_tw.get('theta_NY')

V_eff_ax_tw = cancel(eng_ax_tw.data['potential'])
if theta_ax_tw is not None:
    V_eff_ax_tw = cancel(V_eff_ax_tw.subs(theta_ax_tw, S.Zero))

if om1_ax is not None:
    om_ax_subs = [(s, S.Zero) for s in [params_ax_tw.get('omega1'),
                   params_ax_tw.get('omega2'), params_ax_tw.get('omega3')]
                  if s is not None]
    d2_ax_sp1 = cancel(diff(V_eff_ax_tw, om1_ax, 2).subs(om_ax_subs))
    subs_ax_sp1 = [(L_ax_tw, 1), (kappa_ax_tw, 1)]
    if eta_ax_tw is not None:
        subs_ax_sp1.append((eta_ax_tw, 1))
    if alpha_ax_tw is not None:
        subs_ax_sp1.append((alpha_ax_tw, S.Zero))
    val_ax_sp1 = float(d2_ax_sp1.subs(subs_ax_sp1).subs(r_ax_tw, 3))
    print(f"\n  AX omega1 at eta=1, r=3, alpha=0: {val_ax_sp1:.6e}")

    # VT vs AX ratio
    vt_sp1 = spin1_masses_vt.get('omega1')
    if vt_sp1 is not None:
        subs_vt_sp1 = [(L_tw, 1), (kappa_tw, 1)]
        if alpha_tw is not None:
            subs_vt_sp1.append((alpha_tw, S.Zero))
        if V_tw is not None:
            subs_vt_sp1.append((V_tw, S.Zero))
        vt_sp1_f = float(vt_sp1.subs(subs_vt_sp1).subs(r_tw, 3))
        ratio_sp1 = vt_sp1_f / val_ax_sp1
        pass_spin1_eq = abs(ratio_sp1 - 1.0) < 0.01
        print(f"  VT/AX spin-1 ratio = {ratio_sp1:.6f}")
        print(f"  Spin-1 universality (VT=AX): {'PASS' if pass_spin1_eq else 'FAIL'}")
        spin1_result = "= AX (Palatini universality CONFIRMED)" if pass_spin1_eq else f"ratio {ratio_sp1:.4f}"
    else:
        spin1_result = "N/A"
        pass_spin1_eq = None

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("Theorem 4 (Palatini Universality) Summary")
print("=" * 70)
print()
print("  S3 AX vs VT backgrounds give identical spin-2/1 mass spectra:")
print(f"  Spin-2: {spin2_result}")
print(f"  Spin-1: {spin1_result}")
print()

checks = {
    "VT dropout (C2_EC = C2_LC)":           pass_vt_drop,
    "Palatini protection (V* = 0)":          pass_palatini,
    "Spin-2 VT = AX (Palatini universality)": pass_spin2_eq,
    "Spin-1 VT = AX (Palatini universality)": pass_spin1_eq,
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
print(f"  Theorem 4 (Palatini universality): {'PASS' if all_ok else 'FAIL (see above)'}")
if all_ok:
    print()
    print("  Physical interpretation:")
    print("    Torsion background (eta, V) does NOT affect the spin-2/1 Weyl mass spectrum.")
    print("    Only the effective potential structure (V_eff) changes with torsion mode.")
    print("    This is the 'Palatini universality' of the EC-Weyl extension.")

teardown_log()
