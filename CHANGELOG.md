# Changelog

## v4 — 2026-11-02

Changes are relative to v3 (2026-06-20).
The paper is retitled **Curvature Components and Restricted Torsion Variations in Einstein–Cartan Minisuperspace** and reorganized around derivative-level curvature separation and specified homogeneous deformation paths.

- Distinguish LC-Weyl, irreducible EC-Weyl, and Ricci-tracefree EC curvature. All three actions retain the Einstein–Cartan scalar term. The known static mixed-torsion correction belongs to the additional nine-component and four-form sectors and is referenced to revised Paper02.
- Include derivatives of contortion. Circle-dependent axial torsion gives a quadratic-gradient term in the Ricci-tracefree model; the old universal Palatini-protection and nonpropagation claims are withdrawn. The assertion of an unchanged full spin-2 kinetic structure is also withdrawn; the coupled kinetic problem at nonzero quadratic-curvature coupling is not established here. Derive the restricted Euler equations, classify the periodic stationary branches including the exceptional mixed branch, and determine their restricted second variations.
- Separate local Nieh–Yan and Pontryagin densities from periodic integrals. Correct the curvature–coframe contraction to one half; remove the internal-dual Pontryagin proxy, inferred nonzero topological number, and chiral-selection interpretation.
- Use a common volume convention and five independent normalized shape directions. The Nil radial slice point has a nonzero shape gradient and circumference derivative. Replace the old quintet mass classification with path coefficients, and explain the redundant old diagonal coordinates and the gradient term in Hessian transformations.
- Identify the scalar-curvature and Weyl contributions to the Nil cross coefficient; the entire coefficient cannot be attributed to Weyl curvature alone. Separate the torus volume-changing path from volume-preserving shape. With circle period one, its pure-AX cross coefficient is $`24\pi^3Lr\eta^2/\kappa^2`$, replacing the former $`48\pi^4Lr\eta^2/\kappa^2`$.
- Show that constant torsion-axis orientation can affect the additional nine-component sector on anisotropic geometry. Pure static mixing cubics vanish on the stated zero-twist family, while mixed field cubics need not vanish. Classify finite-angle extrema and all degeneracies as eigenspheres of the anisotropy matrix.
- Organize the aligned mixed-torsion term by fluctuation order about a specified background, referencing Paper02's static potential. Remove the former alpha-differentiated “cubic vertex”: the fixed coupling is not a fluctuating field leg. The zero-torsion term is quartic, and nonzero backgrounds can produce quadratic and cubic terms.
- Withdraw the full homogeneous Nil minimum, the $`\lvert\kappa^2\theta\rvert=1`$ Nieh–Yan threshold, unqualified spin-1 tachyon statements, and physical-mass interpretations lacking a kinetic and constraint analysis. The spherical shape coefficients and directional twist coefficients, including the Nil anisotropy, are static derivatives in the stated coordinates. On a specified coupled shape–twist family, isolate the additional nine-component derivative term and prove strict positivity of the common quartic frequency coefficient. The general coupled field and constraint problem remains outside the scope.
- Distinguish compact metrics along the flat Nil paths using an integral harmonic-form invariant, while retaining discrete sign identifications and the nonstationarity of the reference point.
- Organize the Japanese and English manuscripts around the curvature components and variational results. Keep the main revision note in Appendix D and the detailed changes here; Appendix C gives the coordinate transformations and path comparisons.
- Replace the old generic engine and historical proof/figure runners with a self-contained supplement for the revised formulas. `script/paper03ec.py` provides exact derivation, explicit recomputation into a new directory, and saved-data rendering. `paper_results.model` exposes the selectors `LC`, `EC_IRR`, and `EC_TRACEFREE` for the restricted circle density.
- Replace the three old figures and their datasets with three figures showing curvature components, Nil shape paths, and coupled periodic profiles. The figure notebook reads the new saved dataset by default. The former visualization notebook and old `dppu` API are not interfaces to v4 results.

## v3 — 2026-06-20

The title was unchanged from v1. A note was added distinguishing the static twist Hessian, used for the AX/VT comparison and the three-topology table, from the physical spin-1 mode; it cited the then-current Paper III tachyonic, non-ghost spin-1 result. The nine-theorem presentation was otherwise unchanged. v4 withdraws the spin-1 tachyon interpretation, as listed above.

## v2 — 2026-04-25

The title was unchanged from v1. Theorem 6 renamed the Nil "EC false vacuum" to an "EC slice minimum" on the zero-torsion section and added a full homogeneous Hessian criterion: local minimum for $`\lvert\kappa^2\theta\rvert<1`$, marginal at equality, saddle above. The Nil 0+2+1+1 splitting was restated on this slice-minimum branch. v4 withdraws the full homogeneous minimum and the threshold, as listed above.

## v1 — 2026-04-06

First release under the title *Homogeneous Three-Topology Comparison and Mode Dictionary in Einstein–Cartan + Nieh–Yan Theory: Geometric Structure of EC-Weyl Coupling*, presenting nine theorems on S³×S¹, T³×S¹ and Nil³×S¹.

The changes to the results and interpretations of v1–v3 are listed under v4 above and summarized in Appendix D of the revised manuscript. Figures, data and saved outputs of v1–v3 use earlier formulas and should not be combined with the v4 formulas.
