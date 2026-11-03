# Curvature Components and Restricted Torsion Variations in Einstein–Cartan Minisuperspace

First release: v1, 2026-04-06. Previous revisions: v2, 2026-04-25; v3, 2026-06-20.
Current version: **v4**, revised **2026-11-02**. [日本語](README_ja.md) · [Changes](CHANGELOG.md) · [Archival DOI](https://doi.org/10.5281/zenodo.19425147)

The paper compares LC-Weyl, irreducible EC-Weyl, and Ricci-tracefree EC curvature while retaining the Einstein–Cartan scalar term. It derives the effects of circle-dependent torsion and analyzes specified static shape, twist, and torsion-axis paths on S³, T³, and Nil³. The additional nine-component and four-form sectors carry the circle torsion-gradient term. Five independent shape directions resolve the Nil gradient and two exactly flat paths, while rotating the torsion axis exposes the orientation dependence of the nine-component sector. Periodic stationary branches of the Ricci-tracefree action at positive coupling and finite-angle orientation extrema are classified. A coupled shape–twist family separates a nine-component quadratic frequency term from a strictly positive quartic curvature coefficient shared by the three actions. An integral harmonic-form invariant distinguishes compact metrics along the flat Nil paths.

- [English PDF](DPPUv2-paper03ec.pdf)
- Japanese manuscript: [introduction](DPPUv2-paper03ec_sec01.md), [geometry and actions](DPPUv2-paper03ec_sec02.md), [circle torsion](DPPUv2-paper03ec_sec03.md), [static deformations](DPPUv2-paper03ec_sec04.md), [orientation and interactions](DPPUv2-paper03ec_sec05.md), [coupled periodic profiles](DPPUv2-paper03ec_sec06.md), [conclusions](DPPUv2-paper03ec_sec07.md)
- Appendices: [projections](DPPUv2-paper03ec_appA.md), [derivatives](DPPUv2-paper03ec_appB.md), [coordinate maps and paths](DPPUv2-paper03ec_appC.md), [reproduction and revision note](DPPUv2-paper03ec_appD.md), [fixed Nil quotient](DPPUv2-paper03ec_appE.md)
- [Supplement instructions](script/README.md), [saved data](data/README.md), [figure notebook](script/DPPUv2_Paper_Figures_EC.ipynb)

Static identities and potential classifications established in [Paper02](https://github.com/Muacca/DPPUv2-paper02) are reference results. Conventions follow revised [Paper01](https://github.com/Muacca/DPPUv2-paper01), and the scalar-action comparison uses revised [Paper03](https://github.com/Muacca/DPPUv2-paper03). Results hold for the specified Euclidean variation families; extensions to general joint profiles and to the field and constraint equations of the full theory are left for future work.

## Reproduce the results

From the `script` directory, with Python 3.11 or newer:

```bash
python -m pip install -r requirements.txt
python -B paper03ec.py verify
python -B paper03ec.py recompute --output ../recomputed
python -B paper03ec.py render --data ../data/paper03ec_v4 --output ../rendered
```

`verify` derives exact identities without writing output. `recompute` evaluates formulas and creates the fixed datasets; `render` reads saved data and does not run the symbolic calculation. Output directories must be new. Earlier datasets and the old generic engine do not reproduce v4 results.

## Build the PDF

```bash
cd LaTeX
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The result is `LaTeX/main.pdf`. The repository-root PDF is the manuscript copy for reading. Use the manuscript and dataset version labels when selecting repository or archival files.

Author: Muacca. Code license: [MIT](LICENSE.txt).
