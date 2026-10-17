# Blocked items

Moved here after three failed attempts. Fix the item (or split it) and move it
back to BACKLOG.md to retry.
- `hemocytometer-chamber-depth` Expose the existing chamber factor through the CLI so 0.2 mm and custom slides are accurate. Acceptance: `--chamber-depth` (and `--chamber-factor`) change cells/mL by the expected ratio and a test covers 0.1 vs 0.2 mm. Demo: `lab-tools hemocytometer count --live 120 --dead 8 --squares 5 --chamber-depth 0.2`. Files: lab_tools/lab_ops/hemocytometer.py, tests/test_hemocytometer.py
  blocked: 3 attempts failed (2026-10-17)
