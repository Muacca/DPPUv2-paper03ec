# Reproducing Paper03EC v4

The supplement is limited to the manuscript's families. It is self-contained and does not import another paper's source. Use Python 3.11 or later and [requirements.txt](requirements.txt) in this directory. The notebook additionally requires a Jupyter environment.

## Entry points

From this `script` directory:

```bash
python -m pip install -r requirements.txt
python -B paper03ec.py verify
python -B paper03ec.py recompute --output ../recomputed
python -B paper03ec.py render --data ../data/paper03ec_v4 --output ../rendered
```

Output paths must be new. The manifest lists hashes of the saved files. Rendering checks hashes and reads CSV data; it does not derive curvature or regenerate the dataset.

Open [DPPUv2_Paper_Figures_EC.ipynb](DPPUv2_Paper_Figures_EC.ipynb) for the notebook interface. The adjacent `paper03ec_saved_data_viewer.py` exposes `show(name, source=...)`. Choices are `curvature_components`, `nil_shape_paths`, and `coupled_profiles`. The notebook displays saved data by default; its explicit `RECOMPUTE` cell is disabled by default and, when enabled, writes to a new time-stamped folder `recomputed-notebook-YYYYMMDD-HHMMSS` next to `script`.

## Calculations and conventions

| File | Reader-facing role |
|---|---|
| `paper_results/cartan.py` | Four-dimensional Cartan curvature, geometric torsion, W/R2/R3/TF projections. |
| `paper_results/exterior.py` | Independent exterior-form curvature, Cartan/Bianchi identities, characteristic densities. |
| `paper_results/derive.py` | Exact circle, orientation, static twist, shape, and coordinate identities. |
| `paper_results/extensions.py` | Exact stationary-branch, finite-angle, coupled-profile and Nil fixed-quotient coefficients. |
| `paper_results/coupled.py` | Analytic coframe jets, frequency curvature coefficients and quadrature for (J-01). |
| `paper_results/model.py` | Numerical restricted circle invariants and the three action densities. |
| `paper_results/data.py` | Explicit generation and integrity-checked loading of fixed figure datasets. |
| `paper_results/figures.py` | Rendering saved datasets. |

Indices 0, 1, 2 are spatial; index 3 is the circle. Torsion is geometric torsion, with `h=eta/r`. Coordinate `y` has period 1, and `z=L*y` is physical distance. Dot arguments are derivatives with respect to `z`. Spatial structure eigenvalues `n` are constant. Metric signature is Euclidean. The norm sums all four indices without a two-form normalization factor.

`model.circle_density` accepts these selectors. All retain the EC scalar and use the negative integrated-action convention.

| Selector | Quadratic invariant |
|---|---|
| `LC` | Squared LC-Weyl curvature |
| `EC_IRR` | Squared irreducible EC-Weyl curvature |
| `EC_TRACEFREE` | Squared Ricci-tracefree EC curvature, including R2 and R3 |

The default density includes scalar and Nieh–Yan total derivatives. `bulk=True` removes their periodic boundary contribution under the constant-geometry contract. Neither function includes the volume factor `c_M*L*r**3`.

For exact calculations, `cartan.curvature(G, C, deriv)` uses `deriv(direction, expression)`. Omitting it means a static connection. `derive.py` supplies the circle derivative explicitly; substituting profiles into static curvature is not a dynamic calculation.

The old `dppu` package, proof scripts, and visualization notebook are not v4 APIs. The current notebook and CLI reproduce the revised manuscript. Nonstationary shape coefficients are not physical masses.

## Coupled periodic coefficients

`coupled.integrate(samples=512, amplitude=0.5, r=3, L=1, eta=1, V=1)` returns
the mean curvature-invariant polynomials in ascending frequency order. The
degree-four coefficients are assembled from analytic first and second coframe
jets before periodic quadrature; no frequency fit is used. For this module,
primes are derivatives with respect to period-one `y`, and `n` is a positive
integer frequency, distinct from the spatial structure eigenvalues used by
`model.circle_components`. Keep all scalar controls and the amplitude fixed.

`coupled.action_polynomial(mean, model="EC_IRR", alpha=0, kappa=1, r=3, L=1)`
includes the four-volume and returns the integrated action polynomial. Supply
the same `r,L` as in the quadrature. The three model selectors are shared with
the circle-density API. Complete NY vanishes on this coupled family. Use
`numpy.polynomial.polynomial.polyval(n, coefficients)` for integer frequencies.
The `n^4` coefficient controls this fixed-amplitude sequence, not all profiles.
