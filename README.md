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

## cfu — worked example

A 1:10 serial dilution plated at 0.1 mL per plate with duplicate counts:

| Dilution factor | Replicate A | Replicate B |
|---:|---:|---:|
| 0.01 (10⁻²) | 245 | 231 |
| 0.001 (10⁻³) | 42 | 38 |

All four counts are inside the USP/EP countable range of 30–300 CFU. For each plate:

```
CFU/mL = colony count / (volume plated in mL × dilution factor)
245 / (0.1 × 0.01)  = 245,000 CFU/mL
231 / (0.1 × 0.01)  = 231,000 CFU/mL
 42 / (0.1 × 0.001) = 420,000 CFU/mL
 38 / (0.1 × 0.001) = 380,000 CFU/mL
```

The reported mean across the four usable plates is:

```
(245,000 + 231,000 + 420,000 + 380,000) / 4 = 319,000 CFU/mL
```

The between-plate CV is about 29.82% (sample standard deviation ≈ 95,118 CFU/mL).

Run the same calculation with the one-line command:

```bash
lab-tools cfu quick --counts 245,231,42,38 --dilutions 0.01,0.001 --replicates 2 --volume 0.1
```

The report prints `CFU/mL: 3.19e+05` and the same uncertainty and dilution details.

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
