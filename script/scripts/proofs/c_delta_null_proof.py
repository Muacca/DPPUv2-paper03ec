"""
C_delta Null Proof
==================

Theorem 3 (delta-sector null theorem):
  In the formal MIXING sector of the unified engine, the auxiliary cubic
  coefficient

      C_delta = d3V_eff / d(delta0)d(delta1)d(delta2) |_(delta=0)

  vanishes identically for all three topologies (S3, T3, Nil3) and for all
  torsion backgrounds (AX, VT, MX). Consequently dC_delta/dalpha = 0.

Interpretation:
  The EC extension introduces no new homogeneous auxiliary cubic channel in the
  delta-sector. New cubic structure appears only in the physical eta-V channel.

Author: Muacca
Date: 2026-03-31
"""

import os
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_LIB_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", ".."))
_DATA_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "..", "data"))
sys.path.insert(0, _LIB_DIR)
os.makedirs(_DATA_DIR, exist_ok=True)

from sympy import diff, S

from dppu.utils.tee_logger import setup_log, teardown_log
from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

setup_log(__file__, log_dir=_DATA_DIR)

print("=" * 70)
print("C_delta Null Proof  -- Theorem 3")
print("=" * 70)
print()


def analyze_case(topology_type, topo_name, torsion_mode):
    cfg = DOFConfig(
        topology=topology_type,
        torsion_mode=torsion_mode,
        fiber_mode=FiberMode.MIXING,
        enable_squash=False,
        enable_shear=False,
        ny_variant=NyVariant.FULL,
    )
    engine = UnifiedEngine(cfg)
    engine.run()

    params = engine.data["params"]
    potential = engine.data["potential"]

    delta0 = params["delta0"]
    delta1 = params["delta1"]
    delta2 = params["delta2"]
    alpha = params["alpha"]

    c_delta = diff(potential, delta0, delta1, delta2).subs(
        {delta0: S.Zero, delta1: S.Zero, delta2: S.Zero}
    )
    dc_dalpha = diff(c_delta, alpha)

    zero_pass = c_delta == S.Zero or getattr(c_delta, "is_zero", False) is True
    alpha_pass = dc_dalpha == S.Zero or getattr(dc_dalpha, "is_zero", False) is True

    print(f"{'-' * 65}")
    print(f"  {topo_name}  [{torsion_mode.value}]")
    print(f"{'-' * 65}")
    print(f"  C_delta = {c_delta}")
    print(f"  dC_delta/dalpha = {dc_dalpha}")
    print(f"  Null check: {'PASS' if zero_pass else 'FAIL'}")
    print(f"  alpha-independence: {'PASS' if alpha_pass else 'FAIL'}")
    print()

    return {
        "topology": topo_name,
        "mode": torsion_mode.value,
        "c_delta": c_delta,
        "dc_dalpha": dc_dalpha,
        "passed": zero_pass and alpha_pass,
    }


results = []
for topo_type, topo_name in [
    (TopologyType.S3, "S3 x S1"),
    (TopologyType.T3, "T3 x S1"),
    (TopologyType.NIL3, "Nil3 x S1"),
]:
    for mode in [Mode.AX, Mode.VT, Mode.MX]:
        results.append(analyze_case(topo_type, topo_name, mode))

print("=" * 70)
print("C_delta Summary  -- Theorem 3")
print("=" * 70)
print()

for topo_name in ["S3 x S1", "T3 x S1", "Nil3 x S1"]:
    print(f"  {topo_name}:")
    for mode_name in ["AX", "VT", "MX"]:
        record = next(r for r in results if r["topology"] == topo_name and r["mode"] == mode_name)
        label = "PASS" if record["passed"] else "FAIL"
        print(
            f"    {mode_name:<3}: C_delta = {record['c_delta']}, "
            f"dC_delta/dalpha = {record['dc_dalpha']}  [{label}]"
        )
    print()

all_pass = all(r["passed"] for r in results)

print("  Physical interpretation:")
print("  The formal auxiliary MIXING sector stays null under the EC extension.")
print("  No new homogeneous cubic channel appears in delta0, delta1, delta2.")
print("  The new EC cubic structure is confined to the physical eta-V channel.")
print()
print(f"  Theorem 3 verification: {'PASS' if all_pass else 'FAIL'}")

teardown_log()
