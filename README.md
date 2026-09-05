# lab-tools

CLI tools for pharmaceutical, microbiology and molecular-biology laboratory
data — the calculations a QC lab does every day, packaged as a single tested
Python library with one entry point.

```bash
pip install .
lab-tools cfu --help
```

## Tools

| Domain | Tool | What it does | Status |
|---|---|---|---|
| Microbiology / QC | `cfu` | CFU enumeration from dilution plates (USP ⟨61⟩/⟨62⟩, EP 2.6.12) | ✅ |
| | `mpn` | Most probable number | planned |
| | `media-fill` | Aseptic process simulation (USP ⟨1116⟩) | planned |
| | `growth-curve` | Logistic / Gompertz growth kinetics | planned |
| | `bioburden-spc` | I-MR control charts (Western Electric) | planned |
| Sterility / thermal | `d-z-f0` | D-value, Z-value, F0 lethality | planned |
| | `sterility` | Sterility-test sample sizes (USP ⟨71⟩) | planned |
| | `endotoxin` | Endotoxin limits + MVD (USP ⟨85⟩) | planned |
| Molecular biology | `qpcr` | ΔΔCt relative gene expression | planned |
| | `primers` | Annealing temperature + dimer analysis | planned |
| | `restriction` | Restriction digest simulation | planned |
| | `phylogeny` | Neighbor-joining trees | planned |
| Analytical | `elisa` | 4PL standard curves | planned |
| | `lod-loq` | LOD / LOQ (ICH Q2(R1)) | planned |
| | `hplc-sst` | HPLC system suitability | planned |
| Lab operations | `hemocytometer` | Cell count / viability | planned |
| | `pipette-cal` | Pipette calibration (ISO 8655) | planned |

## Development

```bash
pip install -e . pytest
python -m pytest -q
```

The project grows through a daily backlog-driven automation (see
`BACKLOG.md` and `scripts/`): one item is implemented and committed per day,
only if the test suite stays green.

## License

MIT © Collins Amatu Gorgerat
