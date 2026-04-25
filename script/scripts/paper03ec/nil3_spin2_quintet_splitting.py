"""
Nil3 Spin-2 Quintet Splitting (Theorem 7)
==========================================

Computes the Nil3 x S1 spin-2 quintet splitting pattern under EC-Weyl.
The 5-plet of spin-2 perturbations (eps, s, q3, q4, q5) splits as:
  LC (isotropic): 5-plet (degenerate, like S3 SO(3) symmetry)
  EC slice minimum: 0 + 2 + 1 + 1 (split by C^2_01 = 1/R asymmetry)

Specifically:
  q3 = 0 (zero mode): hyperbolic rotation in (e0, e1) plane preserves C^2_01
  q4 = q5 (2-plet): e0 and e1 equivalent in C^2_01
  eps != s: squash and shear have different masses

squash-shear 2x2 mass matrix:
  M = [[d2V/deps2, d2V/deps_ds],
       [d2V/deps_ds, d2V/ds2]]

at the Nil3 EC slice minimum (r0 = (4*kappa/sqrt(3))*sqrt(|alpha|), eta=0).

Note: off-diagonal DOFs (q3, q4, q5) require enable_offdiag_shear=True,
which may not be implemented. The script handles both cases gracefully.

Author: Muacca
Date: 2026-03-30
"""

import sys
import os
import math

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_LIB_DIR    = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", ".."))
_DATA_DIR   = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "..", "data"))
sys.path.insert(0, _LIB_DIR)
os.makedirs(_DATA_DIR, exist_ok=True)

from dppu.utils.tee_logger import setup_log, teardown_log
setup_log(__file__, log_dir=_DATA_DIR)

from sympy import cancel, diff, Matrix, Rational, S, sqrt, symbols
import sympy as sp

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.topology.unified import DOFConfig as _DOFConfig
from dataclasses import fields
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

print("=" * 70)
print("Nil3 Spin-2 Quintet Splitting  -- Theorem 7")
print("=" * 70)
print()

# ── Step 1: Nil3 shear implementation check ───────────────────────────────────
print("-" * 70)
print("Step 1: Nil3 squash + shear engine")
print("-" * 70)

cfg_shear = DOFConfig(
    topology=TopologyType.NIL3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    enable_shear=True,
    ny_variant=NyVariant.FULL,
)
eng_shear = UnifiedEngine(cfg_shear)
eng_shear.run()

params_sh = eng_shear.data['params']
R_sh      = params_sh.get('R') or params_sh['r']
eps_sh    = params_sh['epsilon']
s_sh      = params_sh['s']
eta_sh    = params_sh.get('eta')
L_sh      = params_sh['L']
kappa_sh  = params_sh['kappa']
alpha_sym = params_sh.get('alpha')

print(f"  epsilon = {eps_sh}")
print(f"  s       = {s_sh}")

# C^2_01 dependence check
C_array = eng_shear.data.get('structure_constants')
if C_array is not None:
    c201 = C_array[2, 0, 1]
    expected = -(1 + eps_sh)**Rational(-4, 3) * (1 + s_sh)**(-2) / R_sh
    diff_c = cancel(c201 - expected)
    print(f"\n  C^2_01 = {c201}")
    print(f"  C^2_01 matches expected -(1+eps)^(-4/3)*(1+s)^(-2)/R: {diff_c.is_zero is True}")
else:
    print("  (structure_constants not directly accessible)")

V_eff_shear = eng_shear.data['potential']

# ── Step 2: 2x2 squash-shear mass matrix ─────────────────────────────────────
print()
print("-" * 70)
print("Step 2: Spin-2 squash-shear 2x2 mass matrix at Nil3 EC slice minimum")
print("-" * 70)
print()

# Compute second derivatives at eps=s=0
V_eff_s = cancel(V_eff_shear)
d2V_deps2  = cancel(diff(V_eff_s, eps_sh, 2))
d2V_ds2    = cancel(diff(V_eff_s, s_sh, 2))
d2V_depds  = cancel(diff(diff(V_eff_s, eps_sh), s_sh))

d2V_deps2_iso = cancel(d2V_deps2.subs([(eps_sh, 0), (s_sh, 0)]))
d2V_ds2_iso   = cancel(d2V_ds2.subs([(eps_sh, 0), (s_sh, 0)]))
d2V_depds_iso = cancel(d2V_depds.subs([(eps_sh, 0), (s_sh, 0)]))

print(f"  d2V/deps2 |_(eps=s=0) = {d2V_deps2_iso}")
print(f"  d2V/ds2   |_(eps=s=0) = {d2V_ds2_iso}")
print(f"  d2V/deps_ds|_(eps=s=0) = {d2V_depds_iso}")
print()

M_spin2 = Matrix([
    [d2V_deps2_iso, d2V_depds_iso],
    [d2V_depds_iso, d2V_ds2_iso],
])
print(f"  2x2 mass matrix M:")
print(f"    [[{d2V_deps2_iso}, {d2V_depds_iso}],")
print(f"     [{d2V_depds_iso}, {d2V_ds2_iso}]]")

# Eigenvalues (symbolic)
det_M = cancel(M_spin2.det())
tr_M  = cancel(M_spin2.trace())
print(f"\n  det(M) = {det_M}")
print(f"  tr(M)  = {tr_M}")

# Numerical evaluation at EC slice minimum: r0 = 4*kappa/sqrt(3), a=kappa=L=1, eta=0
r0_num    = float(4 / sqrt(3))
kappa_num = 1.0
L_num     = 1.0
alpha_num = -1.0   # a=1

subs_num = [(L_sh, L_num), (kappa_sh, kappa_num), (R_sh, r0_num)]
if eta_sh is not None:
    subs_num.append((eta_sh, 0))
if alpha_sym is not None:
    subs_num.append((alpha_sym, alpha_num))

M_eps2_num  = float(d2V_deps2_iso.subs(subs_num))
M_s2_num    = float(d2V_ds2_iso.subs(subs_num))
M_epds_num  = float(d2V_depds_iso.subs(subs_num))

print()
print(f"  Numerical at r0={r0_num:.4f}, a=kappa=L=1, eta=0 (alpha=-1):")
print(f"    d2V/deps2  = {M_eps2_num:.6e}")
print(f"    d2V/ds2    = {M_s2_num:.6e}")
print(f"    d2V/deps_ds = {M_epds_num:.6e}")

tr_num  = M_eps2_num + M_s2_num
det_num = M_eps2_num * M_s2_num - M_epds_num**2
disc_num = tr_num**2 - 4 * det_num

print()
print(f"  Eigenvalues (numerical):")
if disc_num >= 0:
    lam_p = (tr_num + math.sqrt(disc_num)) / 2
    lam_m = (tr_num - math.sqrt(disc_num)) / 2
    print(f"    lambda_+ = {lam_p:.6e}")
    print(f"    lambda_- = {lam_m:.6e}")
    if lam_p > 0 and lam_m > 0:
        print("    Both eigenvalues > 0: EC slice minimum STABLE in squash+shear")
        mass_stable = True
    elif lam_m < 0:
        print("    WARNING: shear direction has tachyon (lambda_- < 0)")
        mass_stable = False
    else:
        print(f"    lambda_+ = {lam_p:.4e}, lambda_- = {lam_m:.4e}")
        mass_stable = lam_m >= 0
else:
    print(f"    Complex eigenvalues (disc = {disc_num:.4e})")
    mass_stable = False

# Compare squash-only and squash+shear for eps
cfg_sq_only = DOFConfig(
    topology=TopologyType.NIL3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    enable_shear=False,
    ny_variant=NyVariant.FULL,
)
eng_sq = UnifiedEngine(cfg_sq_only)
eng_sq.run()
params_sq = eng_sq.data['params']
R_sq    = params_sq.get('R') or params_sq['r']
eps_sq  = params_sq['epsilon']
V_sq    = eng_sq.data['potential']
d2V_sq  = cancel(diff(V_sq, eps_sq, 2).subs(eps_sq, 0))

subs_sq = [(params_sq['L'], L_num), (params_sq['kappa'], kappa_num), (R_sq, r0_num)]
if params_sq.get('eta') is not None:
    subs_sq.append((params_sq['eta'], 0))
if params_sq.get('alpha') is not None:
    subs_sq.append((params_sq['alpha'], alpha_num))

M_sq_num = float(d2V_sq.subs(subs_sq))
print()
print(f"  Cross-check: squash-only vs squash+shear (eps component):")
print(f"    squash-only d2V/deps2 = {M_sq_num:.6e}")
print(f"    squash+shear d2V/deps2 = {M_eps2_num:.6e}")
match_eps = abs(M_eps2_num - M_sq_num) < abs(M_sq_num) * 1e-6
print(f"    Match: {'PASS' if match_eps else 'FAIL (check shear implementation)'}")

# ── Step 3: Off-diagonal DOFs (q3, q4, q5) ───────────────────────────────────
print()
print("-" * 70)
print("Step 3: Off-diagonal DOFs (q3, q4, q5) -- implementation check")
print("-" * 70)
print()

cfg_fields = {f.name: f for f in fields(_DOFConfig)}
has_offdiag = 'enable_offdiag_shear' in cfg_fields
print(f"  enable_offdiag_shear field: {'present' if has_offdiag else 'NOT present'}")

offdiag_masses_num = {}
offdiag_implemented = False

if has_offdiag:
    try:
        cfg_offdiag = DOFConfig(
            topology=TopologyType.NIL3,
            torsion_mode=Mode.AX,
            enable_squash=False,
            enable_shear=False,
            enable_offdiag_shear=True,
            ny_variant=NyVariant.FULL,
        )
        eng_od = UnifiedEngine(cfg_offdiag)
        eng_od.run()

        params_od = eng_od.data['params']
        z3 = params_od.get('z3')
        z4 = params_od.get('z4')
        z5 = params_od.get('z5')

        if z3 is not None and z3 != S.Zero:
            print("  Off-diagonal DOFs implemented in Nil3")
            offdiag_implemented = True

            V_eff_od  = cancel(eng_od.data['potential'])
            R_od      = params_od.get('R') or params_od['r']
            L_od      = params_od['L']
            kappa_od  = params_od['kappa']
            eta_od    = params_od.get('eta')
            alpha_od  = params_od.get('alpha')

            base_subs_od = [(L_od, L_num), (kappa_od, kappa_num), (R_od, r0_num)]
            if eta_od is not None:
                base_subs_od.append((eta_od, 0))
            if alpha_od is not None:
                base_subs_od.append((alpha_od, alpha_num))

            for zi, zname, qname in [(z3, 'z3', 'q3'), (z4, 'z4', 'q4'), (z5, 'z5', 'q5')]:
                if zi is None:
                    continue
                subs_zi = base_subs_od + [
                    (zj, sp.S.One) for zj in [z3, z4, z5] if (zj is not None and zj is not zi)
                ]
                V_zi = cancel(V_eff_od.subs(subs_zi))
                d2_zi = cancel(diff(V_zi, zi, 2).subs(zi, sp.S.One))
                mass_num = float(d2_zi) / 2
                offdiag_masses_num[qname] = mass_num
                print(f"    m^2({qname}) = {mass_num:.6e}")
        else:
            print("  Off-diagonal DOFs return S.Zero (not implemented for Nil3)")
    except Exception as exc:
        print(f"  Error: {exc}")
else:
    print("  enable_offdiag_shear not available -- skipping off-diagonal DOF analysis")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("Nil3 Spin-2 Quintet Splitting Summary (Theorem 7)")
print("=" * 70)
print()
print("  At the Nil3 EC slice minimum (r0 = 4*kappa/sqrt(3)*sqrt(|alpha|), eta=0):")
print()

if disc_num >= 0:
    print(f"  2x2 (eps, s) block:")
    print(f"    lambda_+ (heavy squash/shear) = {lam_p:.6e}")
    print(f"    lambda_- (light squash/shear) = {lam_m:.6e}")
else:
    print("  2x2 block: complex eigenvalues")

if offdiag_implemented and offdiag_masses_num:
    print()
    print("  Off-diagonal block:")
    for qname, mass in offdiag_masses_num.items():
        print(f"    m^2({qname}) = {mass:.6e}")
else:
    print()
    print("  Off-diagonal (q3, q4, q5): not computed (enable_offdiag_shear unavailable)")
    print("  Expected pattern (Theorem 7):")
    print("    q3 = 0 (zero mode): hyperbolic rotation in (e0,e1) preserves C^2_01")
    print("    q4 = q5 (2-plet): e0 and e1 equivalent roles in C^2_01")

print()
print("  Splitting pattern (Theorem 7):")
print("    LC isotropic: 5-plet (approximate degeneracy)")
print("    EC slice minimum: 0 + 2 + 1 + 1 (broken by C^2_01 = 1/R direction)")
print()
print("  Contrast with S3:")
print("    S3: SO(3) exact 5-plet degeneracy")
print("    Nil3: C^2_01 breaks SO(3) -> splits quintet")
print()
print(f"  2x2 mass stability: {'STABLE' if mass_stable else 'UNSTABLE'}")
print(f"  eps-component cross-check: {'PASS' if match_eps else 'FAIL'}")

teardown_log()
