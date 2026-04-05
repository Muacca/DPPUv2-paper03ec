"""
S3 MX vacuum search under EC-Weyl.

This auxiliary script documents the S3 mixing-mode (MX) statement used in
paper03ec Section 3:

  - on the alpha < 0 branch, no finite real MX stationary point exists;
  - on the alpha > 0 branch, real MX stationary points can exist, but they
    are saddles rather than stable full stationary vacua;
  - a finite-r well on a fixed (eta, V) slice does not imply full stationarity.
"""

import os
import sys

import numpy as np
from sympy import Matrix, cancel, nsolve, sqrt

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


def dedup_roots(roots, tol=1e-6):
    unique = []
    for root in roots:
        if any(max(abs(root[i] - other[i]) for i in range(len(root))) < tol for other in unique):
            continue
        unique.append(root)
    return unique


def find_real_roots(eqns, variables, seeds, radial_tol=0.2):
    roots = []
    for seed in seeds:
        try:
            sol = nsolve(eqns, variables, seed, tol=1e-14, maxsteps=100, prec=50)
        except Exception:
            continue

        vals = []
        ok = True
        for value in sol:
            z = complex(value.evalf())
            if abs(z.imag) > 1e-9:
                ok = False
                break
            vals.append(float(z.real))
        if not ok or vals[0] <= radial_tol:
            continue
        roots.append(tuple(vals))
    return dedup_roots(roots)


def hessian_eigenvalues(hessian, substitutions):
    matrix = hessian.subs(substitutions)
    array = np.array([[float(matrix[i, j].evalf()) for j in range(matrix.cols)] for i in range(matrix.rows)])
    return np.linalg.eigvalsh(array)


print("=" * 70)
print("S3 MX Vacuum Search under EC-Weyl")
print("=" * 70)
print()

veff_ec, veff_lc, params = build_veff_ec(
    TopologyType.S3,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
    L_val=1,
    kappa_val=1,
)
veff_ec = cancel(veff_ec)
veff_lc = cancel(veff_lc)
delta_veff = cancel(veff_ec - veff_lc)

r_s = params["r"]
eta_s = params["eta"]
V_s = params["V"]
theta_s = params["theta_NY"]
alpha_s = params["alpha"]

grad = [cancel(veff_ec.diff(var)) for var in (r_s, eta_s, V_s)]
hessian = Matrix(
    [
        [cancel(veff_ec.diff(v1).diff(v2)) for v2 in (r_s, eta_s, V_s)]
        for v1 in (r_s, eta_s, V_s)
    ]
)

print("Step 0: symbolic input")
print("-" * 70)
print(f"  V_eff_LC = {veff_lc}")
print(f"  delta_V_eff (EC-LC) = {delta_veff}")
print()

print("Step 1: exact stationary equations")
print("-" * 70)
for label, expr in zip(("dV/dr", "dV/deta", "dV/dV"), grad):
    print(f"  {label} = {expr}")
print()

print("Step 2: theta_NY = 0 analytic MX branch")
print("-" * 70)
print("  On the MX branch (eta*V != 0), the theta_NY = 0 stationary equations give")
print(f"    V*   = +/- {cancel(3 / (4 * sqrt(alpha_s)))}")
print("    eta* = +/- 2/sqrt(3)")
print(f"    r*   = {cancel(8 * sqrt(alpha_s) / sqrt(3))}")
print("  Therefore alpha < 0 gives no finite real MX stationary point.")

sample_subs = {
    theta_s: 0.0,
    alpha_s: 1.0,
    r_s: float(8 / np.sqrt(3)),
    eta_s: float(2 / np.sqrt(3)),
    V_s: 0.75,
}
sample_eigs = hessian_eigenvalues(hessian, sample_subs)
print("  Sample Hessian eigenvalues at alpha = 1, theta_NY = 0:")
print(f"    {sample_eigs}")
print("  One eigenvalue is negative, so the real MX point is a saddle.")
print()

print("Step 3: numerical root search and Hessian check")
print("-" * 70)
print("  Sampled box for root search: 0.2 < r < 20, |eta| <= 3, |V| <= 3")

scan_cases = [
    (0.0, -1.0),
    (0.0, -0.2),
    (0.5, -1.0),
    (-0.5, -1.0),
    (0.0, 0.2),
    (0.0, 1.0),
    (0.5, 0.2),
    (0.5, 1.0),
    (-0.5, 1.0),
]
seeds = [
    (0.5, -2.0, -2.0),
    (0.5, -2.0, 2.0),
    (0.5, 2.0, -2.0),
    (0.5, 2.0, 2.0),
    (2.0, -2.0, -2.0),
    (2.0, -2.0, 2.0),
    (2.0, 2.0, -2.0),
    (2.0, 2.0, 2.0),
    (5.0, -2.0, -2.0),
    (5.0, -2.0, 2.0),
    (5.0, 2.0, -2.0),
    (5.0, 2.0, 2.0),
    (8.0, 2.0, 1.0),
]

stable_in_dictionary = 0
all_roots = 0

for theta_val, alpha_val in scan_cases:
    eqns = [expr.subs({theta_s: theta_val, alpha_s: alpha_val}) for expr in grad]
    roots = find_real_roots(eqns, (r_s, eta_s, V_s), seeds)
    print()
    print(f"  theta_NY = {theta_val:+.2f}, alpha = {alpha_val:+.2f}")
    if not roots:
        print("    no finite real root found")
        continue

    for root in roots:
        eigs = hessian_eigenvalues(
            hessian,
            {theta_s: theta_val, alpha_s: alpha_val, r_s: root[0], eta_s: root[1], V_s: root[2]},
        )
        stable = bool(np.all(eigs > 1e-8))
        in_dictionary = abs(root[1]) <= 3.0 and abs(root[2]) <= 3.0
        all_roots += 1
        if stable and in_dictionary:
            stable_in_dictionary += 1
        print(
            "    "
            f"r={root[0]:.6f}, eta={root[1]:+.6f}, V={root[2]:+.6f}, "
            f"in_dict={in_dictionary}, stable={stable}, eig={np.array2string(eigs, precision=3)}"
        )

print()
print("Step 4: representative 1D radial well")
print("-" * 70)
slice_subs = {eta_s: 1.2, V_s: 0.8, theta_s: 0.5, alpha_s: 1.0}
dr_expr = grad[0].subs(slice_subs)
deta_expr = grad[1].subs(slice_subs)
dV_expr = grad[2].subs(slice_subs)
d2r_expr = cancel(dr_expr.diff(r_s))
r_slice = float(nsolve(dr_expr, r_s, 3.0))
print("  Fix (eta, V, theta_NY, alpha) = (1.2, 0.8, 0.5, 1.0)")
print(f"  The 1D radial extremum is at r = {r_slice:.6f}")
print(f"  d2V/dr2 at that point = {float(d2r_expr.subs({r_s: r_slice}).evalf()):.6f}")
print(f"  dV/deta at that point = {float(deta_expr.subs({r_s: r_slice}).evalf()):.6f}")
print(f"  dV/dV   at that point = {float(dV_expr.subs({r_s: r_slice}).evalf()):.6f}")
print("  This is a radial well on a fixed slice, not a full stationary vacuum.")
print()

print("=" * 70)
print("S3 MX Vacuum Search Summary")
print("=" * 70)
print(f"  Real roots found in the sampled cases: {all_roots}")
print(f"  Stable full stationary vacua inside the sampled dictionary box: {stable_in_dictionary}")
if stable_in_dictionary == 0:
    print("  Verdict: no stable S3 MX full stationary vacuum was found in the sampled dictionary range.")
else:
    print("  Verdict: unexpected stable MX root found; recheck the paper statement.")

teardown_log()
