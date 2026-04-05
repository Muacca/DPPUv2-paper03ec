"""
T3 MX Vacuum Search under EC-Weyl
===================================

Searches for stable vacua in the T3 x S1 mixing-mode (MX) sector
under EC-Weyl coupling (alpha != 0).

Key results:
  LC (alpha=0): no finite MX critical point
  alpha<0 branch:
    on the theta_NY=0 slice, the stationary solutions
      V* = +/- 3/(4*sqrt(alpha)),   eta* = +/- r/(4*sqrt(alpha))
    are purely imaginary, so no finite real MX critical point exists
    delta_V_eff = -256*pi^4*V^2*alpha*eta^2*r/3 > 0 for alpha<0
    => V_eff_EC > V_eff_LC in the MX sector

Conclusion: T3 MX has triple triviality under EC-Weyl:
  - AX dropout (C2_EC = C2_LC = 0)
  - VT dropout (C2_EC = C2_LC = 0)
  - no finite real MX critical point for alpha<0

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

import numpy as np
from scipy.optimize import minimize

from sympy import cancel, lambdify, symbols, sqrt

from dppu.topology.unified import TopologyType
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.action.ec_action import build_veff_ec

print("=" * 70)
print("T3 MX Vacuum Search under EC-Weyl")
print("=" * 70)
print()

# ── Step 0: Build V_eff_EC (T3 MX) and lambdify ──────────────────────────────
print("-" * 70)
print("Step 0: Build V_eff_EC (T3 MX)")
print("-" * 70)
print()

print("  [Building T3 MX engine...]", flush=True)
veff_ec, veff_lc, params = build_veff_ec(
    TopologyType.T3,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
    L_val=1,
    kappa_val=1,
)
delta_veff = cancel(veff_ec - veff_lc)
print("  Done")
print()

print(f"  V_eff_LC (T3 MX) = {cancel(veff_lc)}")
print(f"  delta_V_eff (EC-LC) = {delta_veff}")
print()

r_s     = params.get('r',        symbols('r',        positive=True))
eta_s   = params.get('eta',      symbols('eta',      real=True))
V_s     = params.get('V',        symbols('V',        real=True))
theta_s = params.get('theta_NY', symbols('theta_NY', real=True))
alpha_s = params.get('alpha',    symbols('alpha',    real=True))

# Lambdify for numerical minimization (kappa=L=1 already substituted)
veff_ec_num = lambdify([r_s, eta_s, V_s, theta_s, alpha_s],
                        cancel(veff_ec), 'numpy')
veff_lc_num = lambdify([r_s, eta_s, V_s, theta_s, alpha_s],
                        cancel(veff_lc), 'numpy')


# ── Helper functions ──────────────────────────────────────────────────────────

def classify_mode(eta, V, tol=0.1):
    if abs(eta) < tol and abs(V) < tol:
        return "trivial (eta=V=0)"
    if abs(eta) < tol:
        return "VT (eta~0)"
    if abs(V) < tol:
        return "AX (V~0)"
    return "MX (eta!=0, V!=0)"


def f_ec(x, theta, alpha):
    r, eta, V = float(x[0]), float(x[1]), float(x[2])
    if r <= 0.01:
        return 1e10
    return float(veff_ec_num(r, eta, V, theta, alpha))


def f_lc(x, theta, alpha):
    r, eta, V = float(x[0]), float(x[1]), float(x[2])
    if r <= 0:
        return 1e10
    return float(veff_lc_num(r, eta, V, theta, alpha))


def dedup_points(pts, tol=0.2):
    unique = []
    for pt in pts:
        is_dup = False
        for upt in unique:
            if (abs(pt[0] - upt[0]) < tol and abs(pt[1] - upt[1]) < tol
                    and abs(pt[2] - upt[2]) < tol):
                is_dup = True
                break
        if not is_dup:
            unique.append(pt)
    return unique


# ── Step 1: LC reference (alpha=0) ───────────────────────────────────────────
print("-" * 70)
print("Step 1: LC reference (alpha=0)")
print("-" * 70)
print()
print("  Analytical: V_eff_LC(T3,MX) = (16*pi^4/kappa^2) * r * [V^2*r^2/3 + 2*V*eta*r*theta + 3*eta^2]")
print("  => no finite minimum: dV/deta=0 at eta=0, dV/dr=0 needs V=0 => trivial only")
print()

# ── Step 2: Analytic MX stationary conditions on the theta_NY=0 slice ────────
print("-" * 70)
print("Step 2: Analytic MX stationary conditions on the theta_NY=0 slice")
print("-" * 70)
print()
V_star = 3 * sqrt(1 / alpha_s) / 4
eta_star = r_s * sqrt(1 / alpha_s) / 4
print("  From dV/deta = dV/dV = 0 at theta_NY = 0:")
print(f"    V* = ±{cancel(V_star)}")
print(f"    eta* = ±{cancel(eta_star)}")
print("  Therefore alpha < 0 gives purely imaginary stationary points.")
print("  => no finite real MX critical point exists on the alpha<0 branch")
print()

# ── Step 3: Auxiliary numerical scan on the alpha<0 branch ───────────────────
print("-" * 70)
print("Step 3: Auxiliary numerical scan on the alpha<0 branch")
print("-" * 70)
print()
print("  delta_V_eff = -256*pi^4*V^2*alpha*eta^2*r/3")
print("  alpha<0: delta_V_eff > 0 => V_eff_EC > V_eff_LC in the MX sector")
print()

scan_params = [
    (0.0,  -0.05), (0.0,  -0.1), (0.0,  -0.2),
    (0.5,  -0.05), (0.5,  -0.1), (0.5,  -0.2),
    (1.0,  -0.05), (1.0,  -0.1), (1.0,  -0.2),
    (-1.0, -0.05), (-1.0, -0.1), (-1.0, -0.2),
]

ec_results = []
new_mx_count = 0

for theta, alpha in scan_params:
    local_pts = []
    for r0 in [0.5, 1.0, 2.0, 3.0]:
        for eta0 in [-0.5, 0.5, 1.0, -1.0]:
            for V0 in [-0.5, 0.5, 1.0, -1.0]:
                try:
                    res = minimize(f_ec, [r0, eta0, V0], args=(theta, alpha),
                                   method='Nelder-Mead',
                                   options={'xatol': 1e-8, 'fatol': 1e-10,
                                            'maxiter': 15000})
                    r_opt, eta_opt, V_opt = res.x
                    if r_opt > 0.01 and r_opt < 20:
                        mode = classify_mode(eta_opt, V_opt)
                        local_pts.append((r_opt, eta_opt, V_opt, res.fun, mode))
                except Exception:
                    pass

    unique_pts = dedup_points(local_pts)
    best = min(unique_pts, key=lambda x: x[3]) if unique_pts else None

    if best and "MX" in best[4]:
        new_mx_count += 1
    ec_results.append({'theta': theta, 'alpha': alpha, 'best': best})

print(f"  {'theta_NY':>8}  {'alpha':>7}  {'best minimum':>32}  mode")
print("  " + "-" * 66)
for res in ec_results:
    theta, alpha = res['theta'], res['alpha']
    best = res['best']
    if best:
        r, eta, V, E, mode = best
        flag = "  *** NEW MX ***" if "MX" in mode else ""
        print(f"  {theta:>8.2f}  {alpha:>7.3f}  r={r:.2f} eta={eta:.2f} V={V:.2f} E={E:.4f}  [{mode}]{flag}")
    else:
        print(f"  {theta:>8.2f}  {alpha:>7.3f}  (minimization failed)")

print()
print(f"  Nontrivial MX minima found in sampled alpha<0 configurations: {new_mx_count}/{len(scan_params)}")
print()

# ── Step 4: Hessian check for any sampled MX points ──────────────────────────
print("-" * 70)
print("Step 4: Hessian stability check for sampled MX candidates")
print("-" * 70)
print()

mx_results = [r for r in ec_results if r['best'] and "MX" in r['best'][4]]

if mx_results:
    print(f"  MX points found in {len(mx_results)} parameter configurations")
    for res in mx_results[:3]:
        theta, alpha = res['theta'], res['alpha']
        r_opt, eta_opt, V_opt, E_opt, mode = res['best']
        print(f"\n  theta={theta:.2f}, alpha={alpha:.3f}:")
        print(f"    minimum: r={r_opt:.4f}, eta={eta_opt:.4f}, V={V_opt:.4f}, E={E_opt:.6f}")

        # Numerical Hessian
        h = 1e-4
        x0 = np.array([r_opt, eta_opt, V_opt])
        f_local = lambda x: f_ec(x, theta, alpha)
        H = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                x1 = x0.copy(); x1[i] += h; x1[j] += h
                x2 = x0.copy(); x2[i] += h; x2[j] -= h
                x3 = x0.copy(); x3[i] -= h; x3[j] += h
                x4 = x0.copy(); x4[i] -= h; x4[j] -= h
                H[i, j] = (f_local(x1) - f_local(x2) - f_local(x3) + f_local(x4)) / (4 * h**2)
        eigvals = np.linalg.eigvalsh(H)
        print(f"    Hessian eigenvalues: {eigvals}")
        stable = all(ev > 0 for ev in eigvals)
        print(f"    Stable (all eigenvalues > 0): {'YES' if stable else 'NO'}")
else:
    print("  No nontrivial MX minima found in the sampled alpha<0 branch")
    print()
    print("  Analytical explanation:")
    print(f"    delta_V_eff (EC - LC) = {delta_veff}")
    print("    alpha < 0 => delta_V > 0 everywhere in MX sector")
    print("    => V_eff_EC > V_eff_LC, consistent with the absence of real critical points")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print("=" * 70)
print("T3 MX Vacuum Search Summary")
print("=" * 70)
print()

if new_mx_count == 0:
    conclusion = "T3 MX has no finite real critical point on the alpha<0 branch"
    verdict = "CONFIRMED: T3 is EC-trivial (no new phases under EC-Weyl)"
else:
    conclusion = f"T3 MX has {new_mx_count} stable point(s) — requires further investigation"
    verdict = "UNEXPECTED: further analysis needed"

print(f"  Auxiliary scan: nontrivial MX minima in alpha<0 sample = {new_mx_count}/{len(scan_params)}")
print(f"  Conclusion: {conclusion}")
print()
print(f"  Verdict: {verdict}")
print()
print("  paper03ec Section 4.3 content:")
if new_mx_count == 0:
    print("    On the alpha<0 branch, the theta_NY=0 stationary solutions are not real,")
    print("    and the EC-Weyl correction raises V_eff in the MX sector:")
    print("    delta_V = +256|alpha|V^2*eta^2*r/3 > 0.")
    print("    Combined with AX/VT dropout, T3 exhibits triple triviality.")

teardown_log()
