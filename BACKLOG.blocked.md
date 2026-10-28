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
