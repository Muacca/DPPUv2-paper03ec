# Saved data for Paper03EC v4

`paper03ec_v4/` supplies the three manuscript figures.

| File | Contents |
|---|---|
| `exact_results.json` | Exact projected invariants, scalar/characteristic densities, Euler equations, Nil coefficients, and zero residuals of symbolic identities. |
| `circle_profile.csv` | A fixed smooth periodic profile and its additional curvature norms, Nieh–Yan and Pontryagin densities. |
| `nil_paths.csv` | Five normalized exponential shape paths, with `A=B=1` and zero torsion. |
| `nil_coefficients.csv` | First and second path derivatives at the Nil reference point. |
| `extension_results.json` | Stationary-branch and second-variation formulas, coupled metric-derivative identities, and Nil harmonic-lattice formulas. |
| `coupled_coefficients.json` | Scalar controls, quadrature size, and ascending frequency coefficients of the integrated curvature invariants. |
| `coupled_profiles.csv` | Integer-frequency actions at three couplings, and the constant/derivative parts of the tracefree shift. |
| `manifest.json` | Version, profile, shape chart, scope, and hashes of the preceding files. |

The circle example is off shell. The Nil reference point is not stationary in shape. Plotted coefficients are not particle masses. Exact formulas apply to the stated restricted families, not arbitrary joint metric/torsion profiles.

Use [the public CLI](../script/README.md) to render saved data or explicitly recompute into a new folder. Earlier-version logs and figure outputs are not inputs to this dataset.

The current manifest uses `paper03ec.data.v2`. The loader also accepts the earlier v1 data for its two figures; the coupled figure requires the current dataset. The coupled example is off shell and holds amplitude and scalar controls fixed. Its high-frequency growth is not a lower-bound or physical-stability result for the full action.
