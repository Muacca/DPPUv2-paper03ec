"""
Nil3 EC slice-minimum stability proof.

This script refines paper03ec Theorem 6 by separating:

1. the analytic eta=V=0 slice minimum on Nil3 for alpha < 0, and
2. the full homogeneous stability criterion in (R, eta, V).

Main symbolic result:

  det H_(eta,V) = (262144*pi^8*L^2*alpha^2/9) * (1 - (kappa^2*theta_NY)^2)

Hence the slice stationary point is a full local minimum iff
|kappa^2 * theta_NY| < 1, marginal at equality, and a saddle otherwise.
"""

import os
import sys

import numpy as np
from sympy import Matrix, Rational, cancel, diff, factor, pi, simplify, solve

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_LIB_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", ".."))
_DATA_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "..", "data"))
sys.path.insert(0, _LIB_DIR)
os.makedirs(_DATA_DIR, exist_ok=True)

from dppu.utils.tee_logger import setup_log, teardown_log

setup_log(__file__, log_dir=_DATA_DIR)

from dppu.action.ec_action import build_veff_ec
from dppu.topology.unified import TopologyType
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant


print("=" * 70)
print("Nil3 EC Slice-Minimum Stability Proof")
print("=" * 70)
print()

veff_ec, _, params = build_veff_ec(
    TopologyType.NIL3,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
)

R_s = params.get("R", params["r"])
eta_s = params["eta"]
V_s = params["V"]
theta_s = params["theta_NY"]
alpha_s = params["alpha"]
kappa_s = params["kappa"]
L_s = params["L"]

print("Step 1: eta=V=0 slice potential")
print("-" * 70)
slice_veff = cancel(veff_ec.subs({eta_s: 0, V_s: 0, theta_s: 0}))
dV_dR = cancel(diff(slice_veff, R_s))
r0_solutions = solve(dV_dR, R_s)
r0_sym = r0_solutions[0]
V0_sym = simplify(cancel(slice_veff.subs(R_s, r0_sym)))

print(f"  V_eff|_(eta=V=theta=0) = {slice_veff}")
print(f"  dV/dR = {dV_dR}")
print(f"  r0 = {r0_sym}")
print(f"  V0 = {V0_sym}")
print()

print("Step 2: full homogeneous Hessian at the slice stationary point")
print("-" * 70)
vars_h = (R_s, eta_s, V_s)
H = Matrix(
    [
        [cancel(diff(veff_ec, v1, v2)) for v2 in vars_h]
        for v1 in vars_h
    ]
)
H_at_r0 = simplify(H.subs({R_s: r0_sym, eta_s: 0, V_s: 0}))
spin0_block = H_at_r0[1:3, 1:3]
H_rr = simplify(H_at_r0[0, 0])
minor_eta = simplify(spin0_block[0, 0])
det_spin0 = factor(simplify(spin0_block.det()))
expected_det_exact = factor(
    Rational(262144, 9)
    * pi**8
    * L_s**2
    * alpha_s**2
    * (1 - kappa_s**4 * theta_s**2)
)

print(f"  H_rr = {H_rr}")
print(f"  H_(eta,eta) = {minor_eta}")
print(f"  det H_(eta,V) = {det_spin0}")
print()

det_ok = simplify(det_spin0 - expected_det_exact) == 0

print(
    "  Full local minimum criterion:"
    " H_rr > 0, H_(eta,eta) > 0, det H_(eta,V) > 0"
)
print("  => iff |kappa^2 * theta_NY| < 1")
print("     = 1 : marginal")
print("     > 1 : saddle")
print(f"  determinant closed-form match: {'PASS' if det_ok else 'FAIL'}")
print()

print("Step 3: numerical benchmark")
print("-" * 70)
numeric_cases = [
    (0.5, "local minimum"),
    (1.0, "marginal"),
    (1.2, "saddle"),
]

benchmark_ok = True
for theta_val, expected_label in numeric_cases:
    H_num = H_at_r0.subs(
        {
            alpha_s: -1,
            kappa_s: 1,
            L_s: 1,
            theta_s: theta_val,
        }
    ).evalf()
    eigs = sorted(float(ev) for ev in H_num.eigenvals().keys())
    if eigs[0] > 1e-8:
        actual_label = "local minimum"
    elif abs(eigs[0]) <= 1e-8:
        actual_label = "marginal"
    else:
        actual_label = "saddle"
    benchmark_ok &= actual_label == expected_label
    print(
        f"  theta_NY={theta_val:+.2f}: eigs={np.array2string(np.array(eigs), precision=3)}"
        f"  => {actual_label}"
    )

print()
print("=" * 70)
print("Nil3 EC Slice-Minimum Stability Summary")
print("=" * 70)
print(f"  Slice-potential stationary point exists for alpha<0: {'PASS' if bool(r0_solutions) else 'FAIL'}")
print(f"  Full Hessian criterion det H_(eta,V) ~ 1-(kappa^2 theta_NY)^2: {'PASS' if det_ok else 'FAIL'}")
print(f"  Numerical benchmark (min / marginal / saddle): {'PASS' if benchmark_ok else 'FAIL'}")

teardown_log()
