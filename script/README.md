# DPPUv2 Computation Engine v5 — Scripts Directory

⇒ [日本語](README_ja.md)

**Paper**: "Einstein-Cartan + Nieh-Yan Theory on Homogeneous 3-Manifolds: Unified Analysis of S³×S¹, T³×S¹, and Nil³×S¹" (paper03-EC)

Python packages and execution scripts for numerical and symbolic computation in the Einstein-Cartan + Nieh-Yan framework across three homogeneous topologies (S³×S¹, T³×S¹, Nil³×S¹).

---

## Directory Structure

```
script/
├── docs/                      # Technical documentation and conventions
├── dppu/                      # Main Python package (DPPUv2 Engine v5)
│   ├── geometry/              # Metric, volume form, structure constants
│   ├── connection/            # Levi-Civita connection, Contortion, EC connection
│   ├── curvature/             # Riemann, Ricci, Hodge dual, Pontryagin, Weyl
│   ├── torsion/               # Torsion modes, Ansatz, Nieh-Yan density
│   ├── action/                # Lagrangian, effective potential, stability classification
│   ├── topology/              # Unified engine (S³×S¹, T³×S¹, Nil³×S¹)
│   ├── engine/                # Computation pipeline, logging, checkpointing
│   ├── kk/                    # Kaluza-Klein mode extraction and validation
│   └── utils/                 # Common utilities (Levi-Civita symbol, symbolic computation, visualization)
│
└── scripts/                   # Execution scripts
    ├── paper03ec/             # paper03-EC analysis scripts
    ├── proofs/                # Analytic and symbolic proof scripts
    └── visualize/             # Interactive viewer and notebooks
```

### `docs/` — Documentation

Technical documentation and conventions:
- [DPPUv2 Engine CONVENTIONS](docs/CONVENTIONS.md) — Engine core conventions and specifications
- [DPPUv2 SymPy guideline](docs/SymPy_guideline.md) — SymPy usage guidelines and best practices

---

## Package Overview (dppu/)

| Module | Role | Key Classes / Functions |
|--------|------|------------------------|
| [`geometry`](dppu/geometry/README.md) | Metric and frame field definitions | `build_metric`, `frame_field` |
| [`connection`](dppu/connection/README.md) | EC connection construction | `levi_civita`, `contortion`, `ec_connection` |
| [`curvature`](dppu/curvature/README.md) | Curvature tensors, Pontryagin, Weyl | `RiemannTensor`, `compute_pontryagin_inner_product`, `WeylTensor` |
| [`torsion`](dppu/torsion/README.md) | Torsion structure | `Mode`, `NyVariant`, `build_torsion_tensor` |
| [`action`](dppu/action/README.md) | Action and stability analysis | `build_lagrangian`, `classify_stability` |
| [`topology`](dppu/topology/README.md) | Unified engine for all topologies and DOFs | `UnifiedEngine`, `DOFConfig`, `TopologyType`, `FiberMode` |
| [`engine`](dppu/engine/README.md) | 15-step computation pipeline | `BaseFrameEngine`, `ComputationLogger`, `CheckpointManager` |
| [`kk`](dppu/kk/README.md) | Kaluza-Klein mode extraction | `extract_maxwell`, `extract_mass`, `extract_cs` |
| [`utils`](dppu/utils/README.md) | Common utilities | `epsilon_symbol`, `prove_zero`, `set_style` |

---

## Script Overview (scripts/)

### paper03ec/ — paper03-EC Analysis

#### Mode Dictionary and Properties

| Script | Description |
|--------|-------------|
| `nil3_mode_dictionary.py` | Nil³ mode properties: AX/VT dropout, spin-2 quintet splitting, spin-1 mass splitting, MX background |
| `t3_mode_dictionary.py` | T³ mode properties: flat geometry, α-independence, torsion-volume coupling |

#### MX Vacuum Search

| Script | Description |
|--------|-------------|
| `nil3_mx_vacuum_search.py` | MX vacuum search in Nil³×S¹ parameter space |
| `s3_mx_vacuum_search.py` | MX vacuum search in S³×S¹ parameter space |
| `t3_mx_vacuum_search.py` | MX vacuum search in T³×S¹ parameter space |

#### Spectral Analysis and Theorems

| Script | Description |
|--------|-------------|
| `nil3_spin2_quintet_splitting.py` | Nil³ spin-2 quintet splitting (q₃ = zero mode) |
| `s3_vt_spin_masses.py` | Palatini universality check: AX/VT give identical spin-2/1 spectra (V\*=0 protection) |
| `squash_shear_cross_term.py` | ε-s cross-term decomposition across topologies (T³=kinematic, S³=0, Nil³=curvature) |

#### EC Vertex and EFT

| Script | Description |
|--------|-------------|
| `ec_cubic_vertex.py` | EC cubic vertex: θ-cubic preservation and α-cubic (EC-specific) channel across topologies |
| `nil3_ec_slice_minimum_eft.py` | Nil³ EC slice-minimum EFT: α-induced minimum on η=V=0 |

---

### proofs/ — Analytic and Symbolic Proofs

| Script | Theorem / Content |
|--------|------------------|
| `ax_vt_dropout_proof.py` | **Theorem 1 (AX/VT Dropout)**: C²_EC = C²_LC — Weyl-squared is inactive for single-component torsion |
| `c_delta_null_proof.py` | **Theorem 3 (δ-sector null)**: C_δ ≡ 0 for the MIXING sector across all topologies and modes |
| `palatini_protection_proof.py` | **Palatini protection**: V\* = 0 on-shell for VT background (dV/dV = 0) |
| `t3_flatness_null_test.py` | **T³ flatness**: C²_LC ≡ 0 — no isolated EC slice minimum in flat topology |
| `gamma_scaling_proof.py` | **Theorem 6 (γ=1/2 scaling)**: r₀(α) ∝ √\|α\| and V_eff^c(α) ∝ √\|α\| for Nil³ EC slice minimum |
| `nil3_slice_minimum_stability.py` | **Theorem 6 (Nil³ EC slice minimum)**: full homogeneous Hessian criterion at $(r_0,0,0)$, with minimum / marginal / saddle classification by $|\kappa^2\theta_{\mathrm{NY}}|$ |

---

### visualize/ — Interactive Viewer and Notebooks

| File | Contents |
|------|----------|
| `DPPUv2_Paper_Figures_EC.ipynb` | Jupyter notebook for generating publication figures |
| `DPPUv2_visualize_notebook_v4ec.ipynb` | visualization notebook  |
| `DPPUv2_interactive_viewer_v4ec.py` | Interactive Jupyter widget: topology/mode/NY-variant selector, phase diagram, potential V_eff(r), axis scale toggling |

---

## Quick Start

### Install Dependencies

```bash
pip install sympy numpy scipy pandas tqdm matplotlib
```

### Run Engine with UnifiedEngine

```python
from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

# S³×S¹ with squash, AX torsion, full Nieh-Yan:
cfg = DOFConfig(
    topology=TopologyType.S3,
    enable_squash=True,
    torsion_mode=Mode.AX,
    ny_variant=NyVariant.FULL,
)
engine = UnifiedEngine(cfg)
engine.run()

Veff = engine.data['potential']
fp   = engine.get_free_params()   # dict of active SymPy Symbols
```

### Legacy engine preset names

```python
cfg = DOFConfig.from_engine('S3S1Engine')    # S³×S¹ squash
cfg = DOFConfig.from_engine('T3S1Engine')    # T³×S¹
cfg = DOFConfig.from_engine('Nil3S1Engine')  # Nil³×S¹ squash
cfg.torsion_mode = Mode.MX
engine = UnifiedEngine(cfg)
engine.run()
```

### Run Proof Scripts

```bash
# Theorem 1: AX/VT dropout (C²_EC = C²_LC)
python scripts/proofs/ax_vt_dropout_proof.py

# Theorem 3: delta-sector null (C_δ ≡ 0 for MIXING)
python scripts/proofs/c_delta_null_proof.py

# Palatini protection (V* = 0 for VT background)
python scripts/proofs/palatini_protection_proof.py

# T³ flatness (C²_LC ≡ 0, no isolated EC slice minimum)
python scripts/proofs/t3_flatness_null_test.py

# Theorem 6: Nil³ EC slice-minimum γ=1/2 scaling
python scripts/proofs/gamma_scaling_proof.py

# Theorem 6: Nil³ EC slice minimum stability
python scripts/proofs/nil3_slice_minimum_stability.py
```

### Run Analysis Scripts

```bash
# Mode dictionaries
python scripts/paper03ec/nil3_mode_dictionary.py
python scripts/paper03ec/t3_mode_dictionary.py

# Spectral analysis
python scripts/paper03ec/nil3_spin2_quintet_splitting.py
python scripts/paper03ec/s3_vt_spin_masses.py
python scripts/paper03ec/squash_shear_cross_term.py

# MX vacuum search (all topologies)
python scripts/paper03ec/nil3_mx_vacuum_search.py
python scripts/paper03ec/s3_mx_vacuum_search.py
python scripts/paper03ec/t3_mx_vacuum_search.py
```

---

## paper03-EC Main Results and Corresponding Scripts

| Result | Content (Proof Type) | Script |
|--------|----------------------|--------|
| **Theorem 1: AX/VT Dropout** | C²_EC = C²_LC — single-component torsion does not activate Weyl — **symbolic proof** | `proofs/ax_vt_dropout_proof.py` |
| **Theorem 3: δ-sector null** | C_δ ≡ 0 for MIXING sector across all topologies — **symbolic proof** | `proofs/c_delta_null_proof.py` |
| **Palatini protection** | V\* = 0 on-shell for VT background — **analytic** | `proofs/palatini_protection_proof.py`, `paper03ec/s3_vt_spin_masses.py` |
| **T³ flatness** | C²_LC ≡ 0 — no isolated EC slice minimum for flat topology — **numerical + symbolic** | `proofs/t3_flatness_null_test.py`, `paper03ec/t3_mode_dictionary.py` |
| **Nil³ spin-2 quintet splitting** | q₃ = zero mode; 4+1 splitting — **numerical** | `paper03ec/nil3_spin2_quintet_splitting.py`, `paper03ec/nil3_mode_dictionary.py` |
| **Nil³ EC slice-minimum EFT** | α-induced minimum at η=V=0; γ=1/2 scaling — **analytic + numerical** | `paper03ec/nil3_ec_slice_minimum_eft.py`, `proofs/gamma_scaling_proof.py` |
| **Nil³ EC slice-minimum stability** | full homogeneous Hessian at $(r_0,0,0)$; local minimum for $|\kappa^2\theta_{\mathrm{NY}}|<1$, marginal at equality, saddle above — **analytic + symbolic** | `proofs/nil3_slice_minimum_stability.py` |
| **ε-s cross-term decomposition** | T³=kinematic / S³=0 / Nil³=curvature — **symbolic** | `paper03ec/squash_shear_cross_term.py` |
| **EC cubic vertex** | θ-cubic preserved; α-cubic channel EC-specific — **symbolic** | `paper03ec/ec_cubic_vertex.py` |
| **MX vacuum existence** | Vacuum structure in MX mode per topology — **numerical scan** | `paper03ec/nil3_mx_vacuum_search.py`, `paper03ec/s3_mx_vacuum_search.py`, `paper03ec/t3_mx_vacuum_search.py` |

---

## Topology Reference

| Topology | Lie Group | Background Curvature | Koszul Formula |
|----------|-----------|----------------------|----------------|
| S³×S¹ | SU(2) | +24/r² | Simplified (bi-invariant) |
| T³×S¹ | U(1)³ | 0 (flat) | Trivial |
| Nil³×S¹ | Heisenberg group | −1/(2R²) | General (bi-invariant does not hold) |

---

## Checkpoint Feature

Long-running computations can be interrupted and resumed:

```bash
# Run with checkpointing
python scripts/paper03ec/nil3_mx_vacuum_search.py \
    --output-dir output/ --checkpoint-dir ./checkpoints

# Resume after interruption (simply re-specify the checkpoint directory)
python scripts/paper03ec/nil3_mx_vacuum_search.py \
    --output-dir output/ --checkpoint-dir ./checkpoints
```

---

## License

See LICENSE file in the repository root.

---

**Author**: Muacca
**Version**: DPPUv2 Engine v5
**Date**: 2026-03
