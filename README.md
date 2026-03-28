# Ab initio derivation of the Berry-Keating correction coefficient from the Conrey-Snaith triple correlation formula

**Author:** David Alarcon, Universidad Pablo de Olavide, Sevilla, Spain

**Zenodo DOI:** [10.5281/zenodo.19268985](https://doi.org/10.5281/zenodo.19268985)

## Abstract

The gap ratio statistic of Riemann zeta zeros converges to the GUE prediction
as <r>(T) = R_inf + c/log^2(T) with c_emp = 1.245 +/- 0.040. We decompose c
into two channels using the Conrey-Snaith formula for the triple correlation
R_3 and the conditional gap probability E^cond from the Schur complement of
the sine kernel. The delta R_3 channel, computed ab initio via Richardson
extrapolation at L = 1000, 3000, 10000, and 30000, contributes c_R3 = -1.6,
while the gap probability channel contributes c_E = +2.8 (inferred by
difference), yielding c = c_R3 + c_E = +1.245. Separately, we construct
the Berry-Keating kernel and show that its Fredholm determinant gives the
wrong sign, establishing that the kernel phase is the fundamental obstruction.
The Conrey-Snaith route bypasses this because R_3 = |K|^2 is phase-immune.

## Repository structure

```
main.tex              Main manuscript (10 pages)
references.bib        Bibliography (7 references)
cover_letter.tex      Cover letter
figures/
  fig1_b_rich_map         Richardson coefficient b(v1,v2) heatmap
  fig2_dR3L_convergence   dR3*L convergence (L=1000 vs 3000)
  fig4_b_rich_by_pair     Mean b for different Richardson pairs
  fig6_fredholm_dsigma    Fredholm K^BK: wrong sign diagnostic
  fig7_methods_comparison Comparison of all methods
data/
  dataset_v6_21pts.dat              Empirical dataset (21 points)
  e4_fine_grid_with_Econd.json      Fine grid (1082 pts) with E_cond
  e4_grilla_b_v2_results.json       CS grid (436 pairs, L=18-34)
  e4_anclas_multi_L_results.json    Multi-L anchors (29 pts)
  e4_anclas_v3_multi_L_results.json V3 anchors (29 pts)
src/
  01_conrey_snaith_R3.py          Section 2: R3_GUE, expansion
  02_richardson_extrapolation.py  Section 3: fine grid, heatmap, convergence
  03_two_channel_decomposition.py Section 4: c_R3 with E_cond, sensitivity
  04_fredholm_diagnostic.py       Section 5: sigma_GUE vs sigma_BK, sign
  05_obstructions_analysis.py     Section 6: 5 layers of obstruction
  e4_grilla_b_v2.py               CS grid computation
  e4_anclas_multi_L.py            Multi-L anchor computation
  e4_c_parcial_desde_anclas.py    Integration from anchors
  platt_zeros.py                  Platt zero reader
```

## Compile

```bash
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

## License

CC BY 4.0
