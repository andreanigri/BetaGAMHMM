# Real-data analyses

Annual-state Beta-GAM HMM applications to HMD female death shares by age.

Official cases:

1. Russia 1960–2014, ages 0–40;
2. Italy 1960–2014, ages 0–40;
3. Italy 1950–2019, ages 0–40 (extended-period sensitivity analysis).

The HMM time unit is year. For year `t`, the emission is the 41-dimensional age profile. Ages are not flattened into pseudo-time.

The primary state-number rule selects `lambda` by minimum EDF-based AIC within each `K`, then selects `K` by minimum EDF-based BIC. ICL is retained as a secondary criterion. Occupancy and precision diagnostics do not filter state-number selection.

Robustness analyses include basis dimension, relaxed precision bounds, and a same-data denominator-aware model.

Final inference separates:

- selection-stability bootstrap with re-selection of `K` and `lambda`;
- conditional parametric bootstrap at the selected model for curves, precision parameters, and transition probabilities.

Conditional bootstrap failures receive one lightweight retry with 8 starts instead of 6; successful first-pass replicates are retained unchanged.

Data are not distributed. See `data/README.md` and `../docs/DATA.md`.
