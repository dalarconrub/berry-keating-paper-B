# Ab initio derivation of the Berry-Keating correction coefficient

**Paper:** "Ab initio derivation of the Berry-Keating correction coefficient
from the Conrey-Snaith triple correlation formula"

**Author:** David Alarcon, Universidad Pablo de Olavide, Sevilla, Spain

**Submitted to:** Journal of Number Theory

## Main result

Two-channel decomposition: c = c_R3 + c_E = -1.25 + 2.50 = +1.245

First derivation of the Berry-Keating correction coefficient from first
principles, using the Conrey-Snaith formula for R_3(v_1, v_2; L) with
Richardson extrapolation at L = 1000, 3000, 10000, and 30000.

## Repository structure

```
├── main.tex / main.pdf       Manuscript (7 pages)
├── references.bib             16 references
├── figures/                   7 figures (PDF + PNG)
├── data/
│   ├── dataset_v6_21pts.dat           Empirical dataset (21 points)
│   ├── e4_grilla_b_v2_results.json    CS grid (436 pairs, L=18-34)
│   ├── e4_anclas_multi_L_results.json Multi-L anchors (29 pts, L=15-80)
│   └── e4_anclas_v3_multi_L_results.json  V3 anchors (29 pts, L=20-100)
├── notebooks/                 Analysis notebooks
└── src/                       Computation scripts
```

## Compile

```bash
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

## Related papers

- **Paper A**: Gap ratio statistics of Riemann zeros — Comm. Math. Phys. — [GitHub](https://github.com/dalarconrub/berry-keating-paper-A)
- **Paper B** (this paper): Ab initio derivation of the Berry-Keating correction coefficient — J. Number Theory
- **Paper C**: Berry-Keating spectral convergence rates and the Riemann Hypothesis — Annals of Mathematics — [GitHub](https://github.com/dalarconrub/berry-keating-paper-C)
- **Paper D**: Empirical proof that Berry-Keating convergence implies the Riemann Hypothesis — Nature — [GitHub](https://github.com/dalarconrub/berry-keating-paper-D)
- **Paper E**: Spectral gap functions bounded below by band-limited functions — J. Fourier Anal. Appl. — [GitHub](https://github.com/dalarconrub/berry-keating-paper-E)
- **Data & code**: [GitHub](https://github.com/dalarconrub/berry-keating-riemann)
