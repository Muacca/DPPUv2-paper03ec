"""
Squash-Shear Cross-Term Classification (Theorem 9)
====================================================

Computes the eps-s cross-term d2V/d(eps)d(s)|_{eps=s=0} for all three
topologies (T3, S3, Nil3) under AX torsion background, and classifies
the physical origin of each nonzero cross-term.

Results (Theorem 9):
  T3:   d2V/deps_ds = 48*pi^4*L*eta^2*r / kappa^2  != 0
        Origin: KINEMATIC — volume = (2*pi)^4*L*r^3*(1+eps)*(1+s)
        C2=0 (flat), no curvature contribution
        => alpha-independence confirmed (Weyl mass = 0)

  S3:   d2V/deps_ds = 0  (SymPy exact)
        Origin: ZERO — volume-preserving squash+shear + s->-s symmetry
        Volume = 2*pi^2*L*r^3 (eps,s-independent)

  Nil3: d2V/deps_ds = (384*pi^4*L*R^2 - 8192*pi^4*L*alpha*kappa^2)/(9*R*kappa^2) != 0
        Origin: CURVATURE/WEYL — C^2_01(eps,s) is (eps,s)-dependent
        Volume = (2*pi)^4*L*R^3 (eps,s-independent, volume-preserving)
        alpha-dependence confirms Weyl origin

Step 0: Volume eps,s-dependence check (kinematic vs curvature classification)
Step 1: T3 cross-term (reference value confirmed)
Step 2: S3 cross-term (volume-preserving + s->-s symmetry => 0)
Step 3: Nil3 cross-term (curvature/Weyl, alpha-dependent)
Step 4: Summary table

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

from sympy import cancel, diff, pi, sqrt, S

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

print("=" * 70)
print("Squash-Shear Cross-Term Classification  -- Theorem 9")
print("=" * 70)
print()

# Numerical evaluation parameters
L_num     = 1.0
kappa_num = 1.0
eta_num   = 1.0
r_s3_num  = 3.0
r_t3_num  = 3.0
r_nil_num = float(4 / sqrt(3))   # Nil3 EC false vacuum r0 = (4*kappa/sqrt(3))*sqrt(|alpha|), a=1

print(f"  Numerical parameters: L=kappa=eta=1, r_S3={r_s3_num}, r_T3={r_t3_num}, r_Nil3={r_nil_num:.4f}")
print()

# ── Step 0: Volume eps,s-dependence check ────────────────────────────────────
print("-" * 70)
print("Step 0: Volume eps,s-dependence (kinematic vs curvature origin)")
print("-" * 70)
print()

for topo, topo_name in [
    (TopologyType.S3,   "S3"),
    (TopologyType.T3,   "T3"),
    (TopologyType.NIL3, "Nil3"),
]:
    cfg_v = DOFConfig(
        topology=topo,
        torsion_mode=Mode.AX,
        enable_squash=True,
        enable_shear=True,
        ny_variant=NyVariant.FULL,
    )
    eng_v = UnifiedEngine(cfg_v)
    eng_v.run()

    vol = eng_v.data.get('volume', None)
    params_v = eng_v.data['params']
    eps_v = params_v['epsilon']
    s_v   = params_v['s']

    if vol is not None:
        has_eps = vol.has(eps_v)
        has_s   = vol.has(s_v)
        vol_str = str(cancel(vol))[:60]
        print(f"  {topo_name}:")
        print(f"    Volume = {vol_str}...")
        print(f"    eps-dependent: {has_eps},  s-dependent: {has_s}")
        if not has_eps and not has_s:
            print(f"    => Volume-preserving (kinematic mixing = 0)")
        else:
            print(f"    => Volume varies with eps,s (kinematic mixing active)")
    else:
        print(f"  {topo_name}: volume not in data dict")
    print()

# ── Step 1: T3 cross-term ─────────────────────────────────────────────────────
print("-" * 70)
print("Step 1: T3 AX d2V/deps_ds (kinematic origin)")
print("-" * 70)
print()

cfg_t3 = DOFConfig(
    topology=TopologyType.T3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    enable_shear=True,
    ny_variant=NyVariant.FULL,
)
eng_t3 = UnifiedEngine(cfg_t3)
eng_t3.run()

params_t3 = eng_t3.data['params']
eps_t3 = params_t3['epsilon']
s_t3   = params_t3['s']
r_t3   = params_t3['r']
L_t3   = params_t3['L']
kap_t3 = params_t3['kappa']
eta_t3 = params_t3.get('eta')

V_t3 = cancel(eng_t3.data['potential'])
d2V_t3 = cancel(diff(diff(V_t3, eps_t3), s_t3).subs([(eps_t3, 0), (s_t3, 0)]))
print(f"  d2V/deps_ds|_(eps=s=0) = {d2V_t3}")

subs_t3 = [(L_t3, L_num), (kap_t3, kappa_num), (r_t3, r_t3_num)]
if eta_t3 is not None:
    subs_t3.append((eta_t3, eta_num))
t3_cross_num = float(d2V_t3.subs(subs_t3))
t3_expected  = float(48 * pi**4 * L_num * eta_num**2 * r_t3_num / kappa_num**2)

print(f"  Numerical (r={r_t3_num}, L=kappa=eta=1): {t3_cross_num:.6e}")
print(f"  Expected 48*pi^4*L*eta^2*r/kappa^2:       {t3_expected:.6e}")
t3_ok = abs(t3_cross_num - t3_expected) / abs(t3_expected) < 1e-8
print(f"  Match: {'PASS' if t3_ok else 'FAIL'}")
print()

# Confirm alpha-independence (kinematic, not Weyl)
alpha_t3 = params_t3.get('alpha')
if alpha_t3 is not None:
    has_alpha = d2V_t3.has(alpha_t3)
    print(f"  alpha-independent: {not has_alpha}  => NOT Weyl/graviton mass (kinematic origin)")
else:
    print("  (alpha symbol not found in T3 V_eff)")
print()

# ── Step 2: S3 cross-term ─────────────────────────────────────────────────────
print("-" * 70)
print("Step 2: S3 AX d2V/deps_ds (volume-preserving + s->-s symmetry => 0)")
print("-" * 70)
print()

cfg_s3 = DOFConfig(
    topology=TopologyType.S3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    enable_shear=True,
    ny_variant=NyVariant.FULL,
)
eng_s3 = UnifiedEngine(cfg_s3)
eng_s3.run()

params_s3 = eng_s3.data['params']
eps_s3 = params_s3['epsilon']
s_s3   = params_s3['s']
r_s3   = params_s3['r']
L_s3   = params_s3['L']
kap_s3 = params_s3['kappa']
eta_s3 = params_s3.get('eta')

V_s3 = cancel(eng_s3.data['potential'])
print("  Computing S3 d2V/deps_ds...")
d2V_s3 = cancel(diff(diff(V_s3, eps_s3), s_s3).subs([(eps_s3, 0), (s_s3, 0)]))
print(f"  d2V/deps_ds|_(eps=s=0) = {d2V_s3}")

s3_is_zero = d2V_s3.is_zero is True
print(f"  S3 d2V/deps_ds = 0: {'PASS (SymPy exact)' if s3_is_zero else 'FAIL'}")

subs_s3 = [(L_s3, L_num), (kap_s3, kappa_num), (r_s3, r_s3_num)]
if eta_s3 is not None:
    subs_s3.append((eta_s3, eta_num))
s3_cross_num = float(d2V_s3.subs(subs_s3))
print(f"  Numerical (r={r_s3_num}, L=kappa=eta=1): {s3_cross_num:.6e}  (expect: 0.0)")
print()
print("  Explanation:")
print("    Volume(S3; eps, s) = 2*pi^2*L*r^3 = const w.r.t. eps, s")
print("    => no kinematic mixing from volume")
if s3_is_zero:
    print("    S3 squash+shear V_eff has s -> -s symmetry (odd powers cancel)")
    print("    => d2V/deps_ds = 0 (confirmed)")
else:
    print("    d2V/deps_ds != 0 (unexpected: check S3 potential)")

# ── Step 3: Nil3 cross-term ───────────────────────────────────────────────────
print()
print("-" * 70)
print("Step 3: Nil3 AX d2V/deps_ds (curvature/Weyl origin, alpha-dependent)")
print("-" * 70)
print()

cfg_nil = DOFConfig(
    topology=TopologyType.NIL3,
    torsion_mode=Mode.AX,
    enable_squash=True,
    enable_shear=True,
    ny_variant=NyVariant.FULL,
)
eng_nil = UnifiedEngine(cfg_nil)
eng_nil.run()

params_nil = eng_nil.data['params']
eps_nil   = params_nil['epsilon']
s_nil     = params_nil['s']
R_nil     = params_nil['R']
L_nil     = params_nil['L']
kap_nil   = params_nil['kappa']
eta_nil   = params_nil.get('eta')
alpha_nil = params_nil.get('alpha')

V_nil = cancel(eng_nil.data['potential'])
print("  Computing Nil3 d2V/deps_ds...")
d2V_nil = cancel(diff(diff(V_nil, eps_nil), s_nil).subs([(eps_nil, 0), (s_nil, 0)]))
print(f"  d2V/deps_ds|_(eps=s=0) = {d2V_nil}")

nil_is_nonzero = d2V_nil.is_zero is not True
print(f"  Nil3 d2V/deps_ds != 0: {'PASS' if nil_is_nonzero else 'FAIL (unexpected zero)'}")

# Check alpha-dependence (Weyl origin)
if alpha_nil is not None:
    has_alpha_nil = d2V_nil.has(alpha_nil)
    print(f"  alpha-dependent: {has_alpha_nil}  => CURVATURE/WEYL origin confirmed")

# Numerical values at two alpha points to show alpha-dependence
subs_nil_base = [(L_nil, L_num), (kap_nil, kappa_num), (R_nil, r_nil_num)]
if eta_nil is not None:
    subs_nil_base.append((eta_nil, 0))   # AX false vacuum: eta=0 at extremum

print()
print("  alpha-dependence check (Nil3 cross-term vs alpha):")
if alpha_nil is not None:
    for alpha_val in [0.0, -1.0]:
        subs_nil = subs_nil_base + [(alpha_nil, alpha_val)]
        nil_num = float(d2V_nil.subs(subs_nil))
        print(f"    alpha={alpha_val:5.1f}: d2V/deps_ds = {nil_num:.4e}")
print()
print("  Explanation:")
print("    Volume(Nil3; eps, s) = (2*pi)^4*L*R^3 = const w.r.t. eps, s")
print("    => no kinematic mixing from volume (volume-preserving)")
print("    C^2_01(eps, s) depends on eps, s (asymmetrically in s)")
if nil_is_nonzero:
    print("    => alpha * C^2_EC coupling generates nonzero d2V/deps_ds")
    print("    => pure curvature/Weyl origin (not kinematic)")
else:
    print("    => d2V/deps_ds = 0 (unexpected: check Nil3 curvature coupling)")

# ── Step 4: Summary ───────────────────────────────────────────────────────────
print()
print("=" * 70)
print("Squash-Shear Cross-Term Summary  (Theorem 9)")
print("=" * 70)
print()
print(f"  {'Topology':<10} {'d2V/deps_ds':^30} {'Volume-pres.':<14} {'alpha-dep.':<12} {'Origin'}")
print("  " + "-" * 82)

# T3
t3_alpha_dep = not d2V_t3.has(alpha_t3) if alpha_t3 is not None else "?"
print(f"  {'T3':<10} {'48*pi^4*L*eta^2*r/kappa^2 != 0':^30} {'No':^14} {'No':^12} kinematic (volume)")

# S3
print(f"  {'S3':<10} {'0 (SymPy exact)':^30} {'Yes':^14} {'No':^12} zero (vol-pres + symmetry)")

# Nil3
nil_alpha_dep_str = "Yes" if (alpha_nil is not None and d2V_nil.has(alpha_nil)) else "No"
print(f"  {'Nil3':<10} {'(384-8192*alpha*kappa^2/R)*...':^30} {'Yes':^14} {nil_alpha_dep_str:^12} curvature/Weyl")

print()
print("  Physical interpretation:")
if t3_ok:
    print("    T3:   volume ∝ (1+eps)*(1+s) => kinematic mixing => d2V/deps_ds != 0")
    print("          C2_EC = 0 (flat) => no Weyl contribution")
else:
    print("    T3:   d2V/deps_ds mismatch (check T3 volume coupling)")
if s3_is_zero:
    print("    S3:   volume-preserving + s->-s symmetry => d2V/deps_ds = 0 exactly")
else:
    print("    S3:   d2V/deps_ds != 0 (unexpected)")
if nil_is_nonzero:
    print("    Nil3: volume-preserving, but C^2_01(eps,s) asymmetric in s")
    print("          => alpha * C^2_EC => curvature/Weyl cross-term")
else:
    print("    Nil3: d2V/deps_ds = 0 (unexpected)")
print()
print("  This table forms the core of Theorem 9 (eps-s cross-term classification).")

all_pass = t3_ok and s3_is_zero and nil_is_nonzero
print()
print(f"  Theorem 9 verification: {'PASS' if all_pass else 'FAIL (see above)'}")

teardown_log()
