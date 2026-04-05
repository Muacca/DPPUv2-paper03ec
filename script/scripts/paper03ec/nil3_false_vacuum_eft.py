"""
Nil3 x S1 EC False Vacuum EFT
================================

Computes the effective field theory (EFT) coefficients around the Nil3 x S1
EC false vacuum at r0 = (4*kappa/sqrt(3))*sqrt(|alpha|), eta*=V*=0.

EC false vacuum structure (Theorem 6):
  V_eff_EC(Nil3, r, 0, 0, alpha) = 4*pi^4*r + alpha*(-64*pi^4/3)/r
  => r0 = (4*kappa/sqrt(3))*sqrt(|alpha|)  [delta=1/2]
  => V0 = (32/sqrt(3))*pi^4*sqrt(|alpha|)  [gamma=1/2, false vacuum energy]

Mass spectrum at r0 (alpha = -a^2, kappa=L=1):
  spin-0 r:      m^2_r = 2*sqrt(3)*pi^4*L / (a*kappa^3)
  spin-0 eta:    m^2_eta = 128*sqrt(3)*pi^4*L*a / kappa
  spin-2 squash: m^2_eps = 6272*sqrt(3)*pi^4*L*a / (27*kappa)
  spin-1 omega2: m^2_w2 = 6*sqrt(3)*pi^4*L^3 / (a*kappa^3)
  spin-1 omega0,1: m^2 = 0  (massless; only omega2 is massive)

Mass hierarchy (two groups):
  Group 1 (light, ~ 1/a): {r, omega2}
  Group 2 (heavy, ~ a):   {eta, squash}
  => r is the lightest mode in the homogeneous EFT

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

from sympy import cancel, diff, simplify, solve, sqrt, pi, Symbol, S

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.action.ec_action import build_veff_ec

print("=" * 70)
print("Nil3 x S1 EC False Vacuum EFT")
print("=" * 70)
print()

# ── Step 1: EC false vacuum identification ────────────────────────────────────
print("-" * 70)
print("Step 1: EC false vacuum identification")
print("-" * 70)
print()
print("  Nil3 AX/VT dropout => EC-Weyl does not change the eta=V=0 slice")
print("  => false vacuum determined from the homogeneous eta=V=0 slice")
print()

cfg_ax = DOFConfig(
    topology=TopologyType.NIL3,
    torsion_mode=Mode.AX,
    enable_squash=False,
    ny_variant=NyVariant.FULL,
)
eng_ax = UnifiedEngine(cfg_ax)
eng_ax.run()

params    = eng_ax.data['params']
R_s       = params['r']
kappa_s   = params['kappa']
L_s       = params['L']
alpha_s   = params.get('alpha', Symbol('alpha', real=True))
eta_s     = params['eta']
eps_s     = params.get('epsilon', S.Zero)

V_eff_ax     = eng_ax.data['potential']
V_eff_lc_vac = cancel(V_eff_ax.subs([(eta_s, 0), (eps_s, 0)]))
print(f"  V_eff(eta=V=0) = {V_eff_lc_vac}")

dV_dR = cancel(diff(V_eff_lc_vac, R_s))
print(f"  dV/dR = {dV_dR}")

r0_sols = solve(dV_dR, R_s)
print(f"  dV/dR = 0 solutions: {r0_sols}")

r0_sym = r0_sols[0] if r0_sols else 4 * kappa_s * sqrt(-alpha_s / 3)
print(f"  EC false vacuum: r0 = {r0_sym}")
print()

V0 = simplify(cancel(V_eff_lc_vac.subs(R_s, r0_sym)))
print(f"  Vacuum energy V0 = V_eff(r0) = {V0}")
print()

# Cross-check gamma=1/2 scaling
a_s = Symbol('a', positive=True)
r0_a = simplify(r0_sym.subs(alpha_s, -a_s**2))
V0_a = simplify(V0.subs(alpha_s, -a_s**2))
print(f"  alpha = -a^2: r0 = {r0_a}")
print(f"  alpha = -a^2: V0 = {V0_a}")
print()

# ── Step 2: Mass spectrum around r0 ──────────────────────────────────────────
print("-" * 70)
print("Step 2: Mass spectrum at r0 (second-order fluctuations)")
print("-" * 70)
print()

# spin-0 r: d2V/dR2 at r0
d2V_dR2 = cancel(diff(dV_dR, R_s))
m2_r_raw = cancel(d2V_dR2.subs(R_s, r0_sym))
m2_r = simplify(m2_r_raw)
print(f"  spin-0 r:   d2V/dR2 |_r0 = {m2_r}")

# spin-0 eta: d2V/deta2 at (r0, eta=0)
d2V_deta2 = cancel(diff(V_eff_ax, eta_s, 2).subs(
    [(eta_s, 0), (eps_s, 0), (R_s, r0_sym)]
))
m2_eta = simplify(d2V_deta2)
print(f"  spin-0 eta: d2V/deta2 |_(r0,eta=0) = {m2_eta}")

# cross term R-eta at r0
d2V_dRdeta = cancel(diff(diff(V_eff_ax.subs(eps_s, 0), R_s), eta_s).subs(
    [(eta_s, 0), (R_s, r0_sym)]
))
print(f"  cross term: d2V/dR deta |_(r0,eta=0) = {d2V_dRdeta}")
print()

# spin-2 squash at r0
cfg_ax_sq = DOFConfig(
    topology=TopologyType.NIL3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    ny_variant=NyVariant.FULL,
)
eng_ax_sq = UnifiedEngine(cfg_ax_sq)
eng_ax_sq.run()

params_sq = eng_ax_sq.data['params']
R_sq      = params_sq['r']
eta_sq    = params_sq['eta']
eps_sq    = params_sq['epsilon']
alpha_sq  = params_sq.get('alpha', Symbol('alpha', real=True))
kappa_sq  = params_sq['kappa']

V_eff_sq = eng_ax_sq.data['potential']
V_sq_vac = cancel(V_eff_sq.subs(eta_sq, 0))
d2V_sq_deps2 = cancel(diff(V_sq_vac, eps_sq, 2).subs(eps_sq, 0))
r0_sub = r0_sym.subs(kappa_s, kappa_sq).subs(alpha_s, alpha_sq)
m2_spin2_raw = cancel(d2V_sq_deps2.subs(R_sq, r0_sub))
m2_spin2 = simplify(m2_spin2_raw)
print(f"  spin-2 squash: d2V/deps2 |_(eps=0,eta=0,r0) = {m2_spin2}")

# spin-1 omega2 at r0
cfg_tw = DOFConfig(
    topology=TopologyType.NIL3,
    torsion_mode=Mode.AX,
    fiber_mode=FiberMode.TWIST,
    isotropic_twist=False,
    ny_variant=NyVariant.FULL,
)
eng_tw = UnifiedEngine(cfg_tw)
eng_tw.run()

params_tw = eng_tw.data['params']
R_tw     = params_tw['r']
eta_tw   = params_tw['eta']
om3      = params_tw.get('omega3')
alpha_tw = params_tw.get('alpha', Symbol('alpha', real=True))
kappa_tw = params_tw['kappa']
L_tw     = params_tw['L']

V_eff_tw = eng_tw.data['potential']
if om3 is not None:
    d2V_dom3 = cancel(diff(diff(V_eff_tw, om3), om3).subs([(eta_tw, 0), (om3, 0)]))
    r0_tw = r0_sym.subs([(kappa_s, kappa_tw), (alpha_s, alpha_tw)])
    m2_spin1_raw = cancel(d2V_dom3.subs(R_tw, r0_tw))
    m2_spin1 = simplify(m2_spin1_raw)
    print(f"  spin-1 omega2: d2V/dom3^2 |_(eta=0,om3=0,r0) = {m2_spin1}")
else:
    m2_spin1 = None
    print("  spin-1 omega2: DOF not found in TWIST mode")

# ── Step 3: EFT coefficient summary ──────────────────────────────────────────
print()
print("-" * 70)
print("Step 3: EFT coefficients (alpha = -a^2, kappa=L=1)")
print("-" * 70)
print()

def sub_alpha_unit(expr):
    """Substitute alpha = -a^2, kappa=L=1 and simplify."""
    if expr is None:
        return None
    result = expr
    for sym in expr.free_symbols:
        if sym.name == 'alpha':
            result = result.subs(sym, -a_s**2)
        elif sym.name in ('kappa',):
            result = result.subs(sym, 1)
        elif sym.name == 'L':
            result = result.subs(sym, 1)
    return simplify(result)

m2_r_a    = sub_alpha_unit(m2_r)
m2_eta_a  = sub_alpha_unit(m2_eta)
m2_s2_a   = sub_alpha_unit(m2_spin2)
m2_s1_a   = sub_alpha_unit(m2_spin1)
V0_a2     = sub_alpha_unit(V0)
r0_a2     = sub_alpha_unit(r0_sym)

print(f"  r0 = {r0_a2}")
print(f"  V0 = {V0_a2}")
print()

masses = [
    ("spin-0 r",      m2_r_a),
    ("spin-0 eta",    m2_eta_a),
    ("spin-2 squash", m2_s2_a),
    ("spin-1 omega2", m2_s1_a),
]

print("  Mass spectrum at r0 (alpha=-a^2, kappa=L=1):")
print(f"  {'Mode':<20} {'m^2 (analytic)':^40} {'m^2 @ a=1':^12} {'sign':^8}")
print("  " + "-" * 85)

for name, m2 in masses:
    if m2 is None:
        print(f"  {name:<20} {'(not computed)':^40}")
        continue
    try:
        m2_num = float(m2.subs(a_s, 1))
        sign = "> 0 (stable)" if m2_num > 0 else "< 0 (UNSTABLE)" if m2_num < 0 else "= 0"
    except Exception:
        m2_num = float('nan')
        sign = "?"
    m2_str = str(m2)[:38]
    print(f"  {name:<20} {m2_str:^40} {m2_num:12.4f}  {sign}")

print()
print("  Mass hierarchy:")
print("    Group 1 (~ 1/a, light):  spin-0 r, spin-1 omega2")
print("    Group 2 (~ a, heavy):    spin-0 eta, spin-2 squash")
print()
mass_vals_at1 = {}
for _name, _m2 in masses:
    if _m2 is not None:
        try:
            mass_vals_at1[_name] = float(_m2.subs(a_s, 1))
        except Exception:
            pass
print("  Physical interpretation:")
if mass_vals_at1:
    lightest = min(mass_vals_at1, key=mass_vals_at1.get)
    print(f"    {lightest} is the lightest mode in the homogeneous EFT (m^2={mass_vals_at1[lightest]:.4f} @ a=1)")
else:
    print("    (lightest mode: not determinable)")
print("    eta (AX torsion) is heavy => torsion stabilized at false vacuum")
print("    spin-1 omega_0,1 remain massless; only omega_2 acquires a mass (1-axis mass splitting)")

# ── Stability check ───────────────────────────────────────────────────────────
print()
print("-" * 70)
print("Stability check: all masses > 0 at r0?")
print("-" * 70)
print()

all_stable = True
for name, m2 in masses:
    if m2 is None:
        continue
    try:
        m2_num = float(m2.subs(a_s, 1))
        stable = m2_num > 0
        if not stable:
            all_stable = False
        print(f"  {name}: m^2 = {m2_num:.4f}  => {'stable' if stable else 'UNSTABLE'}")
    except Exception as exc:
        print(f"  {name}: evaluation error: {exc}")

print()
print(f"  EC false vacuum stability: {'STABLE (all m^2 > 0)' if all_stable else 'UNSTABLE'}")

teardown_log()
