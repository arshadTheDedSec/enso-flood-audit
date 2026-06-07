# enso-flood-audit
Reproducible multi-stage audit of ENSO–flood teleconnections
# ENSO–Flood Audit: A Reproducible Multi-Stage Assessment

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20579277.svg)](https://doi.org/10.5281/zenodo.20579277)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the complete, reproducible pipeline and data for:

> Arshad, S. (2026). *A Reproducible Multi-Stage Audit of ENSO–Flood Teleconnections 
> Across 20 Regulated U.S. River Basins: Only Two of Twenty Apparent Signals Survive 
> Rigorous Validation.* Journal of Hydrology.

We test whether ENSO/ONI winter forcing modulates annual flood extremes at 20 major 
U.S. federal dams using a six-stage reproducible audit framework. Of 20 basins tested, 
only 2 survive rigorous validation — demonstrating a 90% attrition rate and showing 
that most apparent ENSO–flood teleconnections disappear under diagnostic scrutiny and 
multiple-testing correction.

---

## Key Results

| Basin | State | n | FDR p | r | Status |
|-------|-------|---|-------|---|--------|
| Dworshak | ID | 59 | 0.000389 | −0.543 | **SURVIVED** |
| Douglas | TN | 75 | 0.030330 | +0.303 | **SURVIVED** |
| Alamo | AZ | 65 | 0.000025 | +0.143 | FRAGILE |
| Detroit | OR | 75 | 0.060133 | −0.262 | MARGINAL |
| Fontana | NC | 71 | 0.064113 | +0.318 | MARGINAL |
| All others (15) | — | — | >0.14 | — | NULL |

FDR correction: Benjamini–Hochberg, α = 0.10, m = 20 basins.

---

## Repository Structure

```
enso-flood-audit/
├── README.md                    ← This file
├── audit_pipeline.py            ← Complete reproducible pipeline (locked)
├── LICENSE                      ← MIT License
├── data/
│   ├── douglas_data.csv         ← Douglas (TN): 75 water-years peak + ONI
│   ├── dworshak_data.csv        ← Dworshak (ID): 59 water-years peak + ONI
│   ├── ONI_FROZEN.csv           ← NOAA CPC monthly ONI 1950–2025 (locked)
│   └── audit_results_v2.csv     ← Full audit results: all 20 basins, FDR-corrected
└── methodology/
    └── AUDIT_PROTOCOL.txt       ← Locked GEV specification and settings
```

---

## Quick Start

### Requirements

```
Python 3.12+
numpy
scipy
pandas
```

Install dependencies:

```bash
pip install numpy scipy pandas
```

### Run the pipeline

```bash
python audit_pipeline.py
```

### Load case study data

```python
import pandas as pd

# Douglas (TN) — positive ENSO coupling
douglas = pd.read_csv('data/douglas_data.csv')

# Dworshak (ID) — negative ENSO coupling
dworshak = pd.read_csv('data/dworshak_data.csv')

# All 20-basin audit results
results = pd.read_csv('data/audit_results_v2.csv')
```

---

## Methodology (Six-Stage Pipeline)

| Stage | Description |
|-------|-------------|
| 1 | Stationary GEV fit (L-BFGS-B, 4 shape initializations) |
| 2 | Non-stationary GEV fit (ONI NDJFM covariate, z-scored) |
| 3 | Likelihood Ratio Test (LRT), D ~ χ²₂, threshold p < 0.05 |
| 4 | Effect-size gate: Pearson \|r\| ≥ 0.15 |
| 5 | Residual diagnostic: Kolmogorov–Smirnov p ≥ 0.05 |
| 6 | Benjamini–Hochberg FDR correction, α = 0.10, m = 20 |

**GEV parameterization (locked):**
- Location: μ(t) = μ₀ + μ₁·ONI_t
- Scale: log σ(t) = σ₀ + σ₁·ONI_t
- Shape ξ: constant, bounded [−0.5, 0.5]
- ONI window: NDJFM (November–March), z-scored
- Random seed: 42

Full specification: see `methodology/AUDIT_PROTOCOL.txt`

---

## Sensitivity Analyses

**Window sensitivity (Table S1):** Pipeline was re-run with DJF and JFM ONI windows.
Both Douglas and Dworshak survive all three windows with nearly identical statistics.

**Effect-size sensitivity (Table S2):** Pipeline was re-run at |r| ≥ 0.10 and |r| ≥ 0.20.
Both Douglas and Dworshak survive all three effect-size thresholds.

---

## Falsification Tests

Three candidate structural explanations were tested and **all three rejected**:

1. **Hydroclimatic thresholds** — basin area, elevation, snow fraction: no correlation with survival status
2. **Operational characteristics** — operator type, storage capacity: no group separation (Twin Basin Paradox)
3. **Flood seasonality** — Mann–Whitney U p = 0.82: NULL basins have higher seasonality than SURVIVED basins

---

## Data Sources

- **Discharge:** USGS National Water Information System (NWIS) RDB export
- **ONI:** NOAA Climate Prediction Center (CPC), monthly 1950–2025
- **Basin attributes:** USGS GAGES-II database (Falcone, 2011)

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{arshad2026enso,
  author  = {Arshad, Sadman},
  title   = {A Reproducible Multi-Stage Audit of {ENSO}--Flood Teleconnections 
             Across 20 Regulated {U.S.} River Basins: Only Two of Twenty 
             Apparent Signals Survive Rigorous Validation},
  journal = {Journal of Hydrology},
  year    = {2026}
}
```

**Data and code deposit:**

```
Arshad, S. (2026). ENSO–Flood Audit: Multi-stage reproducibility testing 
across 20 U.S. federal dams — Data and code. Zenodo. 
https://doi.org/10.5281/zenodo.XXXXXXX
```

---

## License

MIT License. See [LICENSE](LICENSE) file.

---

## Contact

Sadman Arshad  
BUET Civil Engineering  
2304052@ce.buet.ac.bd
