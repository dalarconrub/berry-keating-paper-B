# Ab initio derivation of the Berry-Keating correction coefficient

**Paper:** "Ab initio derivation of the Berry-Keating correction coefficient
from the Conrey-Snaith triple correlation formula"

**Author:** David Alarcon, Universidad Pablo de Olavide, Sevilla, Spain

**Submitted to:** Communications in Mathematical Physics

## Main result

c_teo = 1.23 (99% of c_emp = 1.245 +/- 0.040)

First derivation of the Berry-Keating correction coefficient from first
principles, using the Conrey-Snaith formula for R_3(v_1, v_2; L).

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

## Companion papers

- **Paper 1** (Nature): Empirical measurement + mechanism + (H1)=>RH theorem
- **Paper 2** (Annals): Rigorous proof of the theorem
