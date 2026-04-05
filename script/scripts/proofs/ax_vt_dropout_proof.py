"""
AX/VT Dropout Proof
===================

Theorem 1 (AX/VT dropout): For single-component torsion backgrounds (AX: V=0,
or VT: eta=0), C2_EC = C2_LC for all three topologies. The Weyl-squared
EC coupling is inactive under these backgrounds.

Verification:
  AX mode (V=0): C2_EC = C2_LC  (6 cases: 3 topologies x AX/VT = PASS)
  VT mode (eta=0): C2_EC = C2_LC
  MX mode (V*eta != 0): C2_EC = C2_LC + 16*V^2*eta^2/(3*r^2)

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

from sympy import cancel, diff, simplify, S

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.curvature.ricci import compute_ricci_tensor
from dppu.curvature.weyl import compute_weyl_tensor, compute_weyl_scalar
from dppu.action.ec_action import compute_c2_ec

print("=" * 70)
print("AX/VT Dropout Proof  -- Theorem 1")
print("=" * 70)
print()


def analyze_topology(topology_type, topo_name, torsion_mode, expect_zero=True):
    """Verify AX/VT dropout or MX activation for one topology + mode combination."""
    mode_label = torsion_mode.value
    print(f"\n{'='*65}")
    print(f"  {topo_name}  [{mode_label} mode]")
    print(f"{'='*65}")

    cfg = DOFConfig(
        topology=topology_type,
        torsion_mode=torsion_mode,
        fiber_mode=FiberMode.NONE,
        enable_squash=False,
        enable_shear=False,
        ny_variant=NyVariant.FULL,
    )
    engine = UnifiedEngine(cfg)
    engine.run()
    data   = engine.data
    params = data['params']
    eta    = params['eta']
    alpha  = params['alpha']

    c2_lc = data.get('weyl_scalar', S.Zero)
    c2_ec = compute_c2_ec(data)

    print(f"  C2_LC = {c2_lc}")
    print(f"  C2_EC = {c2_ec}")

    delta_c2 = simplify(c2_ec - c2_lc)
    is_zero = delta_c2 == S.Zero or (hasattr(delta_c2, 'is_zero') and delta_c2.is_zero is True)
    if expect_zero:
        print(f"  C2_EC - C2_LC = {delta_c2}  (expect: 0 for AX/VT)")
        passed = is_zero
        print(f"  Dropout check: {'PASS' if passed else 'FAIL'}")
    else:
        print(f"  C2_EC - C2_LC = {delta_c2}  (expect: nonzero for MX)")
        passed = not is_zero
        print(f"  Activation check: {'PASS' if passed else 'FAIL'}")

    return {
        'topo_name': topo_name,
        'mode':      mode_label,
        'c2_lc':     c2_lc,
        'c2_ec':     c2_ec,
        'delta_c2':  delta_c2,
        'passed':    passed,
    }


# ── AX mode analysis ──────────────────────────────────────────────────────────
print("\n[AX mode: V=0, eta != 0]")
results_ax = {}
for topo_type, topo_name in [
    (TopologyType.S3,   "S3 x S1"),
    (TopologyType.T3,   "T3 x S1"),
    (TopologyType.NIL3, "Nil3 x S1"),
]:
    results_ax[topo_name] = analyze_topology(topo_type, topo_name, Mode.AX, expect_zero=True)

# ── VT mode analysis ──────────────────────────────────────────────────────────
print("\n\n[VT mode: eta=0, V != 0]")
results_vt = {}
for topo_type, topo_name in [
    (TopologyType.S3,   "S3 x S1"),
    (TopologyType.T3,   "T3 x S1"),
    (TopologyType.NIL3, "Nil3 x S1"),
]:
    results_vt[topo_name] = analyze_topology(topo_type, topo_name, Mode.VT, expect_zero=True)

# ── MX mode: verify non-zero delta ────────────────────────────────────────────
print("\n\n[MX mode: V != 0, eta != 0  -- expect C2_EC != C2_LC]")
results_mx = {}
for topo_type, topo_name in [
    (TopologyType.S3,   "S3 x S1"),
    (TopologyType.T3,   "T3 x S1"),
    (TopologyType.NIL3, "Nil3 x S1"),
]:
    results_mx[topo_name] = analyze_topology(topo_type, topo_name, Mode.MX, expect_zero=False)

# ── Summary ──────────────────────────────────────────────────────────────────
topo_names = ["S3 x S1", "T3 x S1", "Nil3 x S1"]

print(f"\n\n{'='*70}")
print("  3-topology summary: AX/VT dropout (Theorem 1)")
print(f"{'='*70}")
print()

for mode_label, results in [("AX", results_ax), ("VT", results_vt)]:
    print(f"  Mode: {mode_label}")
    for name in topo_names:
        r = results[name]
        sym = "PASS" if r['passed'] else "FAIL"
        print(f"    {name:<18}: C2_EC - C2_LC = {r['delta_c2']}  [{sym}]")
    print()

print("  Mode: MX  (expect C2_EC != C2_LC)")
for name in topo_names:
    r = results_mx[name]
    nonzero = not (r['delta_c2'] == S.Zero or
                   (hasattr(r['delta_c2'], 'is_zero') and r['delta_c2'].is_zero is True))
    sym = "CONFIRMED" if nonzero else "UNEXPECTED ZERO"
    print(f"    {name:<18}: delta_C2 = {r['delta_c2']}  [{sym}]")
print()

ax_pass  = all(r['passed'] for r in results_ax.values())
vt_pass  = all(r['passed'] for r in results_vt.values())
mx_nonz  = all(
    not (r['delta_c2'] == S.Zero or
         (hasattr(r['delta_c2'], 'is_zero') and r['delta_c2'].is_zero is True))
    for r in results_mx.values()
)

print(f"  AX dropout (all 3 topologies): {'PASS' if ax_pass else 'FAIL'}")
print(f"  VT dropout (all 3 topologies): {'PASS' if vt_pass else 'FAIL'}")
print(f"  MX non-zero delta (all 3):     {'PASS' if mx_nonz else 'FAIL'}")
print()
print("  Theorem 1 verification: " + ("PASS" if ax_pass and vt_pass else "FAIL"))
print()
print("  Physical interpretation:")
print("  The V x eta cross product in C2_EC = C2_LC + 16*V^2*eta^2/(3*r^2)")
print("  vanishes when either V=0 (AX) or eta=0 (VT).")
print("  Only mixed (MX) backgrounds activate the EC-Weyl coupling.")

teardown_log()
