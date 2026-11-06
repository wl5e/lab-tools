# Blocked items

Moved here after three failed attempts. Fix the item (or split it) and move it
back to BACKLOG.md to retry.
- `qpcr-efficiency-formula` Apply the efficiency-corrected fold-change formula. Acceptance: fold change uses `(1+E)^-ddCq` and equals the current result when every E is 1.0; tested numerically. Demo: `lab-tools qpcr --input qpcr.csv --efficiencies efficiencies.csv`. Files: lab_tools/molecular/qpcr.py, tests/test_qpcr.py
  blocked: 3 attempts failed (2026-10-26)
- `qpcr-efficiency-output` Expose the corrected fold change in the Result and CSV output. Acceptance: the output gains an `efficiency_corrected` column without breaking existing columns; a test asserts the schema. Demo: `lab-tools qpcr --input qpcr.csv --efficiencies efficiencies.csv --format csv`. Files: lab_tools/molecular/qpcr.py, tests/test_qpcr.py
  blocked: 3 attempts failed (2026-10-26)
- `qpcr-grubbs-output` Report flagged/excluded replicates in the qPCR output. Acceptance: the output lists flagged Cq values and the means are recomputed when `--outliers grubbs` is used; tested. Demo: `lab-tools qpcr --input qpcr.csv --outliers grubbs`. Files: lab_tools/molecular/qpcr.py, tests/test_qpcr.py
  blocked: 3 attempts failed (2026-10-27)
- `elisa-4pl-ci` Compute delta-method confidence intervals for 4PL unknowns. Acceptance: `concentration_ci` returns finite bounds for a known curve instead of NaN; tested against a synthetic curve. Demo: `lab-tools elisa unknowns --input samples.csv --ci 95`. Files: lab_tools/analytical/elisa.py, tests/test_elisa.py
  blocked: 3 attempts failed (2026-10-28)
- `elisa-4pl-output` Attach the 4PL confidence interval to unknown-sample output. Acceptance: the text and CSV outputs include lower/upper concentration bounds; a test asserts both. Demo: `lab-tools elisa unknowns --input samples.csv --ci 95 --format csv`. Files: lab_tools/analytical/elisa.py, tests/test_elisa.py
  blocked: 3 attempts failed (2026-10-28)
- `hplc-per-injection-n` Average plate count across replicate injections. Acceptance: per-injection N is computed and the reported N is the mean, not replicate 0; tested on a 3-injection fixture. Demo: `lab-tools hplc-sst --input peaks.csv --injection-col injection`. Files: lab_tools/analytical/hplc_sst.py, tests/test_hplc_sst.py
  blocked: 3 attempts failed (2026-10-28)
- `hplc-sst-signal-to-noise` Add a signal-to-noise column to HPLC SST. Acceptance: a noise input or blank-injection column yields S/N per peak; tested with a known ratio. Demo: `lab-tools hplc-sst --input peaks.csv --noise-col noise`. Files: lab_tools/analytical/hplc_sst.py, tests/test_hplc_sst.py
  blocked: 3 attempts failed (2026-10-30)
- `cfu-ci-output` Report the CFU confidence interval in text and JSON. Acceptance: both outputs include lower/upper CFU/mL; tested. Demo: `lab-tools cfu calculate --counts 30 --dilution 1e-2 --volume 1 --ci`. Files: lab_tools/microbiology/cfu.py, tests/test_cfu.py
  blocked: 3 attempts failed (2026-11-02)
- `cfu-tntc-extrapolation` Add a conservative TNTC lower-bound mode. Acceptance: all-TNTC plates return a documented lower bound instead of an error, and a test covers it. Demo: `lab-tools cfu calculate --counts TNTC --dilution 1e-2 --volume 1 --tntc lower-bound`. Files: lab_tools/microbiology/cfu.py, tests/test_cfu.py
  blocked: 3 attempts failed (2026-11-02)
- `bioburden-cp` Compute process capability Cp against user alert/action limits. Acceptance: `--lsl` and `--usl` produce Cp from the within-subgroup sigma; tested numerically. Demo: `lab-tools bioburden-spc --input data.csv --lsl 0 --usl 100 --capability`. Files: lab_tools/microbiology/bioburden_spc/chart.py, tests/test_bioburden_spc.py
  blocked: 3 attempts failed (2026-11-02)
- `bioburden-cpk` Compute Cpk and report the limiting side. Acceptance: Cpk equals min(CPU, CPL) and the report names the limiting side; tested. Demo: `lab-tools bioburden-spc --input data.csv --lsl 0 --usl 100 --capability`. Files: lab_tools/microbiology/bioburden_spc/chart.py, tests/test_bioburden_spc.py
  blocked: 3 attempts failed (2026-11-03)
- `bioburden-alert-action-limits` Surface alert/action limits on the SPC chart. Acceptance: limits are accepted from CLI or CSV, drawn/labelled, and tested. Demo: `lab-tools bioburden-spc --input data.csv --alert 50 --action 80`. Files: lab_tools/microbiology/bioburden_spc/chart.py, tests/test_bioburden_spc.py
  blocked: 3 attempts failed (2026-11-03)
- `growth-aic` Add SSR and AIC to the growth-curve fit result. Acceptance: AIC is computed from SSR and parameter count, and a test checks a known fit. Demo: `lab-tools growth-curve fit --input od.csv`. Files: lab_tools/microbiology/growth_curve/analyzer.py, tests/test_growth_curve_analyzer.py
  blocked: 3 attempts failed (2026-11-03)
- `growth-bic` Add BIC and per-well model selection. Acceptance: BIC is reported and the selected model is the lower criterion; tested on synthetic logistic/Gompertz data. Demo: `lab-tools growth-curve fit --input od.csv --select bic`. Files: lab_tools/microbiology/growth_curve/analyzer.py, lab_tools/microbiology/growth_curve/main.py, tests/test_growth_curve_analyzer.py
  blocked: 3 attempts failed (2026-11-04)
- `growth-blank-subtraction` Subtract a blank well before fitting. Acceptance: `--blank` subtracts the per-timepoint blank and clamps negatives to zero; tested. Demo: `lab-tools growth-curve fit --input od.csv --blank BLANK`. Files: lab_tools/microbiology/growth_curve/analyzer.py, lab_tools/microbiology/growth_curve/main.py, tests/test_growth_curve_analyzer.py
  blocked: 3 attempts failed (2026-11-04)
- `phylogeny-jc69` Add the Jukes-Cantor distance correction. Acceptance: the corrected distance matches the JC69 formula and rejects p >= 0.75 with a clear error; tested. Demo: `lab-tools phylogeny nj --input distances.phy --distance-model jc69`. Files: lab_tools/molecular/phylogeny.py, tests/test_phylogeny.py
  blocked: 3 attempts failed (2026-11-04)
- `phylogeny-bootstrap` Add bootstrap support values to neighbor-joining. Acceptance: `--bootstrap 100` annotates the Newick tree with support and is deterministic under a fixed seed; tested. Demo: `lab-tools phylogeny nj --input distances.phy --bootstrap 100 --seed 42`. Files: lab_tools/molecular/phylogeny.py, tests/test_phylogeny.py
  blocked: 3 attempts failed (2026-11-06)
- `sterility-membrane-sizing` Add membrane-filtration unit sizing per USP <71>. Acceptance: the membrane mode returns the documented minimum units for a given fill volume and is tested against the compendial table. Demo: `lab-tools sterility size --method membrane --fill-volume 100`. Files: lab_tools/sterility/sterility.py, tests/test_sterility.py
  blocked: 3 attempts failed (2026-11-06)
- `sterility-membrane-cli` Expose `--method {direct,membrane}` on the sterility CLI. Acceptance: both methods run, the default stays direct, and a subprocess test covers the flag. Demo: `lab-tools sterility size --method membrane --fill-volume 100`. Files: lab_tools/sterility/sterility.py, tests/test_sterility.py
  blocked: 3 attempts failed (2026-11-06)
