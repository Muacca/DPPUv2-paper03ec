"""
Nil3 MX vacuum search under EC-Weyl.

This auxiliary script documents the Nil3 mixing-mode (MX) statement used in
paper03ec Section 5:

  - on the MX branch (eta * V != 0), the theta_NY = 0 stationary equations
    have no real solution;
  - on the alpha < 0 branch, numerical root search finds no real MX stationary
    point for sampled nonzero theta_NY values either;
  - the local minimum that remains on alpha < 0 is the eta = V = 0 EC slice minimum,
    not an MX vacuum.
"""

import os
import sys

import numpy as np
from sympy import Matrix, cancel, nsolve, symbols, sqrt

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


def classify_root(root, tol=1e-6):
    if abs(root[1]) < tol and abs(root[2]) < tol:
        return "trivial slice (eta=V=0)"
    if abs(root[1]) < tol:
        return "VT boundary"
    if abs(root[2]) < tol:
        return "AX boundary"
    return "MX"


print("=" * 70)
print("Nil3 MX Vacuum Search under EC-Weyl")
print("=" * 70)
print()

veff_ec, veff_lc, params = build_veff_ec(
    TopologyType.NIL3,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
    L_val=1,
    kappa_val=1,
)
veff_ec = cancel(veff_ec)
veff_lc = cancel(veff_lc)
delta_veff = cancel(veff_ec - veff_lc)

R_s = params["R"]
eta_s = params["eta"]
V_s = params["V"]
theta_s = params["theta_NY"]
alpha_s = params["alpha"]

grad = [cancel(veff_ec.diff(var)) for var in (R_s, eta_s, V_s)]
hessian = Matrix(
    [
        [cancel(veff_ec.diff(v1).diff(v2)) for v2 in (R_s, eta_s, V_s)]
        for v1 in (R_s, eta_s, V_s)
    ]
)

print("Step 0: symbolic input")
print("-" * 70)
print(f"  V_eff_LC = {veff_lc}")
print(f"  delta_V_eff (EC-LC) = {delta_veff}")
print()

print("Step 1: exact stationary equations")
print("-" * 70)
for label, expr in zip(("dV/dR", "dV/deta", "dV/dV"), grad):
    print(f"  {label} = {expr}")
print()

print("Step 2: theta_NY = 0 analytic MX branch")
print("-" * 70)
print("  On the MX branch (eta*V != 0), dV/deta = dV/dV = 0 gives")
print(f"    V*   = +/- {cancel(3 / (4 * sqrt(alpha_s)))}")
print(f"    eta* = +/- {cancel(R_s / (4 * sqrt(alpha_s)))}")
print("  Substituting into dV/dR = 0 yields")
print("    27*R^4 + 12*alpha*R^2 + 64*alpha^2 = 0")
print("  In x = R^2, the discriminant is")
print("    Delta = (12*alpha)^2 - 4*27*64*alpha^2 = -6768*alpha^2 < 0")
print("  Therefore the theta_NY = 0 MX branch has no real stationary point for any alpha != 0.")
print()

print("  For comparison, the alpha < 0 EC slice minimum sits on the trivial slice eta = V = 0.")
alpha_sample = -1.0
R0 = float(4 / np.sqrt(3))
sample_eigs = hessian_eigenvalues(
    hessian,
    {theta_s: 0.0, alpha_s: alpha_sample, R_s: R0, eta_s: 0.0, V_s: 0.0},
)
print(f"    alpha = {alpha_sample:+.1f}, R0 = {R0:.6f}, Hessian eig = {sample_eigs}")
print("  This is the EC slice minimum discussed in Theorem 6, not an MX stationary point.")
print()

print("Step 3: numerical root search for nonzero theta_NY")
print("-" * 70)
print("  Sampled box for root search: 0.2 < R < 20, |eta| <= 3, |V| <= 3")

scan_cases = [
    (0.5, -1.0),
    (-0.5, -1.0),
    (0.5, -0.2),
    (-0.5, -0.2),
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
]

mx_root_count = 0
stable_mx_root_count = 0
trivial_slice_count = 0

for theta_val, alpha_val in scan_cases:
    eqns = [expr.subs({theta_s: theta_val, alpha_s: alpha_val}) for expr in grad]
    roots = find_real_roots(eqns, (R_s, eta_s, V_s), seeds)
    print()
    print(f"  theta_NY = {theta_val:+.2f}, alpha = {alpha_val:+.2f}")
    if not roots:
        print("    no finite real root found")
        continue

    for root in roots:
        label = classify_root(root)
        eigs = hessian_eigenvalues(
            hessian,
            {theta_s: theta_val, alpha_s: alpha_val, R_s: root[0], eta_s: root[1], V_s: root[2]},
        )
        stable = bool(np.all(eigs > 1e-8))
        if label == "MX":
            mx_root_count += 1
            if stable:
                stable_mx_root_count += 1
        if label == "trivial slice (eta=V=0)":
            trivial_slice_count += 1
        print(
            "    "
            f"R={root[0]:.6f}, eta={root[1]:+.6f}, V={root[2]:+.6f}, "
            f"class={label}, stable={stable}, eig={np.array2string(eigs, precision=3)}"
        )

print()
print("=" * 70)
print("Nil3 MX Vacuum Search Summary")
print("=" * 70)
print(f"  MX roots found in sampled nonzero-theta cases: {mx_root_count}")
print(f"  Stable MX roots found in sampled nonzero-theta cases: {stable_mx_root_count}")
print(f"  Trivial-slice roots found in sampled nonzero-theta cases: {trivial_slice_count}")
if stable_mx_root_count == 0:
    print("  Verdict: no stable Nil3 MX stationary point was found; the alpha < 0 minimum remains the eta = V = 0 EC slice minimum.")
else:
    print("  Verdict: unexpected stable MX root found; recheck the paper statement.")

teardown_log()
