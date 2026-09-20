# Bone Marrow Differential

Bone marrow aspirate differential calculator and conservative hematopathology interpretation aid. It computes cell percentages, marrow and non-erythroid blast percentages, the myeloid:erythroid ratio, an approximate age-referenced cellularity description, and pattern-level review flags from entered data.

The project deliberately does **not** treat a differential count as a complete diagnosis. MDS, CMML, plasma-cell neoplasms, aplastic anemia, and genetically defined AML entities require clinical, hematologic, morphologic, cytogenetic, and/or molecular information beyond the fields represented here.

## Main features

- Browser interface using the same Python engine as the CLI via Pyodide; no application server is required.
- Light theme by default with a dark-mode option, responsive layout, compact differential entry, and JSON export.
- CLI for single cases, interactive entry, demonstration cases, and CSV batch processing.
- Strict input checks for negative/non-integer cell counts, invalid percentages, and malformed batch values.
- Separate WHO-HAEM5 and ICC 2022 interpretation notes where terminology or blast thresholds differ.
- Conservative flags for high blast percentages, dysplasia, plasmacytosis, severe hypocellularity, and M:E-ratio shifts without assigning unsupported disease diagnoses.

## Browser application

GitHub Pages deployment is automated from `master` by `.github/workflows/pages.yml`. The live link is added here only after the deployed URL has been verified.

The browser downloads Pyodide and `bone_marrow_differential.py`, then performs analysis locally. Case inputs are not sent to this repository or to an application backend. The interface stores only the selected light/dark theme in browser local storage.

## CLI

Requires Python 3.9 or newer and has no runtime Python dependencies.

```bash
python -m pip install .
bone-marrow-differential --demo normal
bone-marrow-differential --case-id CASE-42 --age 64 --blasts 35 --segs 315 --poly-normo 100 --lymphocytes 50
bone-marrow-differential batch -i sample.csv -o results.csv
```

For JSON output:

```bash
bone-marrow-differential --case-id CASE-42 --age 64 --blasts 35 --segs 315 --poly-normo 100 --lymphocytes 50 --json
```

The CLI does not invent a normal differential when no counts are supplied; missing input is treated as an argument error.

## Python API

```python
from bone_marrow_differential import (
    BoneMarrowCellCounts,
    BoneMarrowDifferentialAnalyzer,
    ClinicalCaseInput,
)

case = ClinicalCaseInput(
    case_id="BM-001",
    patient_age=64,
    counts=BoneMarrowCellCounts(
        blasts=35,
        segmented_neutrophils=315,
        polychromatophilic_normoblasts=100,
        lymphocytes=50,
    ),
    core_cellularity_pct=60,
)

report = BoneMarrowDifferentialAnalyzer.analyze(case)
print(report.to_json())
```

## Interpretation scope

The calculation engine uses entered values as observations, not as proof that diagnostic prerequisites have been met. In particular:

- Dysplasia alone does not establish MDS; cytopenia and exclusion/integrated criteria are required.
- A high marrow plasma-cell percentage does not establish multiple myeloma without clonality and a myeloma-defining event.
- Severe marrow hypocellularity does not establish aplastic anemia without the required peripheral-blood findings and exclusion of alternative causes.
- Monocytosis does not establish CMML without the required peripheral-blood proportion/persistence and supporting or clonal criteria.
- WHO-HAEM5 and ICC 2022 are not collapsed into one AML/MDS threshold rule.
- The displayed IPSS-R blast stratum is reference information only; it is not a complete IPSS-R score.

Outputs are intended for research, education, calculation support, and expert review—not as treatment recommendations or autonomous sign-out.

## Development and testing

```bash
python -m pip install . pytest
python -m pytest -q
python -m compileall -q bone_marrow_differential.py cli.py marrow_mind.py
bone-marrow-differential --demo normal
bone-marrow-differential batch -i sample.csv -o out_smoke.csv
```

Continuous integration tests Python 3.10, 3.11, and 3.12. A browser smoke check loads the page in headless Chrome, waits for Pyodide, executes the Python engine, and verifies a rendered result.

## Technology and browser support

The maintained runtime consists of standard-library Python, HTML, CSS, JavaScript, and Pyodide 314.0.7. Current Chromium-, Firefox-, and WebKit-based browsers with WebAssembly support are expected to work. JavaScript and network access to the pinned Pyodide CDN are required for the browser application; the CLI does not require network access.

## License

MIT License. See [LICENSE](LICENSE).
