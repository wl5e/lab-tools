# Roadmap — October 2026 to January 2027

The daily automation ships one real, tested item per weekday. Weekends are
rest days. `BACKLOG.md` is the executable version of this roadmap; each item
carries its own acceptance criterion, demo command and target files.

## Month 1 — finish what was skipped

The previous daily pipeline never installed `pytest`, so 19 items were checked
off without being implemented. They are re-admitted here, starting with the
atomic ones.

- Expose the hemocytometer chamber factor through the CLI.
- Add JSON output to the bioburden SPC chart.
- Write a real cfu worked example in the README.
- F0: core calculation, time-temperature parser, CLI.
- Primer 3' self-dimer score and CLI.

## Month 2 — quantitative depth and consistency

- qPCR: per-target efficiency, corrected fold change, Grubbs outliers.
- ELISA: 4PL covariance and delta-method confidence intervals.
- HPLC SST: per-injection plates/tailing, resolution, signal-to-noise.
- CFU: Poisson confidence interval, TNTC lower bound.
- Bioburden SPC: Cp/Cpk, alert/action limits.
- Growth curves: AIC/BIC selection, blank subtraction.
- Phylogeny: Jukes-Cantor correction, bootstrap support.
- Sterility: membrane filtration and direct inoculation sizing.
- Endotoxin: IU/unit dosing, EU/unit output, device limits.
- Pipette: evaporation correction, ISO 8655 uncertainty budget.
- LOD/LOQ: signal-to-noise method.
- Media fill: run-size planning.
- Cross-tool: `--json` everywhere, one error format, README index and demos.

## Month 3 — new capabilities and surface polish

- cfu serial-dilution planner and countable-range profiles.
- MPN unequal replicates.
- ELISA 5PL.
- HPLC capacity factor and selectivity.
- F0 target-SAL planning.
- Pipette as-found/as-left modes.
- Route-specific endotoxin K factors.
- Growth-rate comparison.
- Primer GC clamp, circular plasmid digests, phylogeny distance matrices.
- EP 2.6.1 sterility tables.
- Media-fill JSON output.
- A README "What shipped this quarter" section.

## Guardrails

- One item per weekday; weekends off.
- No commit unless the change is real, scoped to at most 3 files, carries a
  test when production code changes, and the full suite is green.
- An item that fails three attempts is moved to `BACKLOG.blocked.md`, not
  silently checked off.
- `RECRUITER-VISIBLE`: CLI exposure, JSON output, README demos and worked
  examples come before internal refactors.
