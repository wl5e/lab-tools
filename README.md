# lab-tools

CLI tools for pharmaceutical, microbiology and molecular-biology laboratory
data — the calculations a QC lab does every day, packaged as a single tested
Python library with one entry point.

```bash
pip install .
lab-tools --help
lab-tools cfu --help
```

## Tools

| Domain | Tool | What it does |
|---|---|---|
| Microbiology / QC | `cfu` | CFU enumeration from dilution plates (USP ⟨61⟩/⟨62⟩, EP 2.6.12) |
| | `mpn` | Most probable number |
| | `media-fill` | Aseptic process simulation (USP ⟨1116⟩) |
| | `growth-curve` | Logistic / Gompertz growth kinetics |
| | `bioburden-spc` | I-MR control charts (Western Electric) |
| Sterility / thermal | `d-z-f0` | D-value, Z-value, F0 lethality |
| | `sterility` | Sterility-test sample sizes (USP ⟨71⟩) |
| | `endotoxin` | Endotoxin limits + MVD (USP ⟨85⟩) |
| Molecular biology | `qpcr` | ΔΔCt relative gene expression |
| | `primers` | Annealing temperature + dimer analysis |
| | `restriction` | Restriction digest simulation |
| | `phylogeny` | Neighbor-joining trees |
| Analytical | `elisa` | 4PL standard curves |
| | `lod-loq` | LOD / LOQ (ICH Q2(R1)) |
| | `hplc-sst` | HPLC system suitability |
| Lab operations | `hemocytometer` | Cell count / viability |
| | `pipette-cal` | Pipette calibration (ISO 8655) |

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
