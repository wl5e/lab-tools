# Blocked items

Moved here after three failed attempts. Fix the item (or split it) and move it
back to BACKLOG.md to retry.
- `qpcr-efficiency-formula` Apply the efficiency-corrected fold-change formula. Acceptance: fold change uses `(1+E)^-ddCq` and equals the current result when every E is 1.0; tested numerically. Demo: `lab-tools qpcr --input qpcr.csv --efficiencies efficiencies.csv`. Files: lab_tools/molecular/qpcr.py, tests/test_qpcr.py
  blocked: 3 attempts failed (2026-10-26)
- `qpcr-efficiency-output` Expose the corrected fold change in the Result and CSV output. Acceptance: the output gains an `efficiency_corrected` column without breaking existing columns; a test asserts the schema. Demo: `lab-tools qpcr --input qpcr.csv --efficiencies efficiencies.csv --format csv`. Files: lab_tools/molecular/qpcr.py, tests/test_qpcr.py
  blocked: 3 attempts failed (2026-10-26)
