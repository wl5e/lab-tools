# Backlog

The daily automation (`scripts/daily.py`) consumes this list top-to-bottom:
it takes the first undone item, implements it with the DeepSeek provider
(`scripts/llm.py`), and ships it only if the full test suite stays green.
Replenishment is handled by `scripts/replenish.py` (propose) and
`scripts/promote.py` (approve).

Keep entries concrete: what changes, and why it matters for lab / GMP work.

## Planned

- [x] `cfu_report_export` Add a formatted human-readable report (not just JSON) to the `cfu` tool. | handler:
- [x] `cfu_edge_cases` Harden `cfu` input validation (zero/negative volumes, empty dilutions) with tests. | handler:
- [x] `cfu_docs` Add a worked example with real numbers to the `cfu` README section. | handler:

- [x] `dzf0-f0-lethality` Add F0 lethality calculation to d-z-f0 — The d-z-f0 tool currently only computes D and Z values, but README advertises F0. Labs validating terminal sterilization need F0 = sum(10^((T-121.1)/Z) * dt) from a time-temperature profile. | handler:
- [x] `primers-self-dimer` Detect 3' primer self-dimers in primers tool — README claims dimer analysis, but primers.py only calculates Tm. A 3' self-complementarity score would catch self-dimer formation, a common cause of PCR failure. | handler:
- [x] `qpcr-efficiency-correction` Support per-target amplification efficiency in qPCR ΔΔCq — Current qpcr hardcodes 2^(-ΔΔCq), assuming 100% efficiency. Adding optional primer efficiency values lets users report corrected fold changes when efficiencies differ from 2. | handler:
- [x] `elisa-4pl-ci` Implement true 4PL inverse concentration confidence intervals — elisa concentration_ci warns and returns NaN for the 4PL model, so unknown concentrations have no uncertainty. Storing the covariance from curve_fit enables a delta-method CI for each sample. | handler:
- [x] `hplc-sst-per-injection` Calculate HPLC SST metrics from every replicate injection — hplc_sst groups peaks but uses only the first row for plate count and tailing, ignoring replicate variation. Per-injection N and tailing should be averaged to give credible system suitability results. | handler:
- [x] `cfu-confidence-interval` Report a 95% confidence interval for CFU estimates — A single CFU/mL value is weak evidence for release or bioburden decisions. Adding a Poisson-based 95% CI gives QC labs an uncertainty window alongside the point estimate. | handler:
- [x] `bioburden-spc-json` Add JSON output mode to bioburden-spc — Bioburden trend data needs to feed LIMS or electronic dashboards; ASCII charts alone are not machine-readable. A --json flag would expose limits, moving ranges, and Western Electric violations for downstream systems. | handler:
- [x] `pipette-cal-evaporation` Apply evaporation correction in pipette calibration — ISO 8655 gravimetric pipette checks are biased when water evaporates from the weighing vessel during measurement. Optional evaporation loss per row should be added back before mass-to-volume conversion. | handler:
- [x] `media-fill-run-size` Calculate required media-fill units for zero-contamination acceptance — Laboratories need to design media fills, not only analyze completed batches. A calculator using binomial/Poisson probability can determine how many units are needed to demonstrate a target contamination rate at a chosen confidence level. | handler:
- [x] `lod-loq-snr-method` Add signal-to-noise based LOD/LOQ calculation — ICH Q2(R1) accepts S/N determination of LOD/LOQ, but lod-loq only uses calibration-curve residual standard deviation. Supporting S/N measurements from blank injections would cover impurity methods where this is the standard approach. | handler:

- [x] `qpcr-outlier-detection` Add Grubbs' outlier detection for qPCR technical replicates — Cq outliers from pipetting or amplification artifacts inflate SD and distort ΔΔCq fold changes. Grubbs' test can flag or exclude them before computing means. | handler:
- [x] `bioburden-spc-capability` Add bioburden process capability (Cp/Cpk) against alert/action limits — SPC charts detect violations but don't quantify how well the process stays within limits. Cpk gives a single number for trend reviews and validation reports. | handler:
- [x] `growth-curve-aic-model-selection` Add AIC/BIC model selection between logistic and Gompertz — R² alone favors the more flexible model; AIC/BIC provides an objective, information-criterion-based choice for each well in growth-curve analysis. | handler:
- [ ] `phylogeny-jukes-cantor` Add Jukes-Cantor distance correction to phylogeny — Neighbor-joining on raw p-distances underestimates evolutionary distances due to multiple substitutions. JC69 correction is a standard first-order fix for distance-based trees. | handler:
- [ ] `sterility-membrane-filtration` Add membrane filtration sample-size calculation to sterility — USP <71> membrane filtration is a common sterility test method with different unit and volume requirements than direct inoculation. Labs need explicit sizing support. | handler:
- [ ] `endotoxin-dose-units` Support unit-based dosing in endotoxin limit calculation — Many biologics are dosed in IU or units, not mg. Converting dose units lets the tool compute correct EU/unit or EU/mL limits per USP <85>. | handler:
- [ ] `pipette-cal-uncertainty-budget` Add ISO 8655 uncertainty budget to pipette calibration — ISO 8655 requires an uncertainty statement alongside pass/fail. A combined and expanded uncertainty tells labs how reliable each volume measurement is. | handler:
- [ ] `hemocytometer-chamber-depth` Support custom hemocytometer chamber depths — The hardcoded 10^4 factor assumes 0.1 mm depth. Chambers of 0.2 mm or custom slides need a corrected factor for accurate cells/mL calculations. | handler:
- [ ] `hplc-sst-signal-to-noise` Add signal-to-noise calculation to HPLC SST — USP <621> requires S/N for some system suitability tests. Adding a noise column or blank injection lets the tool report S/N alongside plates and tailing. | handler:
- [ ] `cfu-tntc-extrapolation` Add conservative TNTC extrapolation option to CFU calculator — When all plates are TNTC, the current tool errors out. Treating TNTC as the countable maximum yields a documented lower-bound estimate instead. | handler:
