# Bone Marrow Differential & Hematopathology Agent

> **Domain:** Hematopathology, Bone Marrow Morphologic Differential, & Myeloid Neoplasm Diagnostics  
> **Standards:** WHO Classification of Haematolymphoid Tumours (5th Edition, 2022) / International Consensus Classification (ICC 2022) / IPSS-R Criteria

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![Tests](https://img.shields.io/badge/Tests-Pytest%2032%20passed-success.svg)
![Standards](https://img.shields.io/badge/Standards-WHO%205th%20%2F%20ICC%202022-brightgreen.svg)

</div>

---

## 📖 Overview

**Bone Marrow Differential Agent** is a specialized clinical hematopathology diagnostic engine and differential analysis pipeline. It evaluates manual 500-cell bone marrow aspirate differentials, calculates myeloid-to-erythroid (M:E) ratios, assesses age-adjusted core biopsy cellularity, quantifies dysplasia across lineages, evaluates Perls' Prussian blue iron stores and ring sideroblasts, and applies WHO 5th Edition (2022) / International Consensus Classification (ICC 2022) criteria for myeloid and plasma cell neoplasms.

---

## 🔬 Hematopathology & Mathematical Formulations

### 1. Manual 500-Cell Aspirate Differential Standard
A manual aspirate smear count of at least **500 nucleated cells** is the international standard to ensure statistical reliability when assessing low-frequency populations such as blasts and plasma cells:

$$\text{Total Nucleated Cells} = \sum_{k} \text{Cell}_k \ge 500$$

### 2. Blast Percentage Calculations
- **Total Nucleated Cell Blast Percentage:**
  $$\text{Blast } \% = \left( \frac{\text{Blasts}}{\text{Total Nucleated Cells}} \right) \times 100$$

- **Non-Erythroid Blast Percentage (FAB Criteria):**
  $$\text{Non-Erythroid Blast } \% = \left( \frac{\text{Blasts}}{\text{Total Nucleated Cells} - \text{Total Erythroid Precursors}} \right) \times 100$$

### 3. Myeloid-to-Erythroid (M:E) Ratio
The ratio of total granulocytic/myeloid lineage precursors to nucleated erythroid precursors:

$$\text{M:E Ratio} = \frac{\text{Blasts} + \text{Promyelocytes} + \text{Myelocytes} + \text{Metamyelocytes} + \text{Bands} + \text{Segs} + \text{Eosinophils} + \text{Basophils}}{\text{Pronormoblasts} + \text{Basophilic Normoblasts} + \text{Polychromatophilic Normoblasts} + \text{Orthochromatophilic Normoblasts}}$$

- **Normal Adult Range:** `1.5:1` to `3.5:1` (typically `2:1` to `4:1`)
- **M:E > 4.5:1:** Myeloid (granulocytic) hyperplasia (infection, leukemoid reaction, CML, G-CSF)
- **M:E < 1.2:1:** Erythroid hyperplasia (hemolysis, hemorrhage recovery, erythropoietin therapy, thalassemia)

### 4. Age-Adjusted Core Cellularity
Expected core biopsy cellularity decreases with age according to the standard clinical regression:

$$\text{Expected Cellularity } (\%) = 100 - \text{Age (years)} \quad (\pm 15\% \text{ normal tolerance window})$$
- *Pediatric adjustments:* Neonates $\approx 95-100\%$, infants/young children $\approx 80-90\%$, minimum floor $\ge 10\%$.

### 5. WHO 5th Edition (2022) / ICC Diagnostic Blast Cutoffs

| Diagnostic Category | Marrow Blast % | Peripheral Blood Blast % | Defining Criteria & Molecular Genetics |
|:---|:---:|:---:|:---|
| **Normal Bone Marrow** | $< 5\%$ | $0\%$ | Normocellular for age, normal M:E ratio ($2:1$ to $4:1$), no dysplasia |
| **MDS with Low Blasts (MDS-LB / MDS-SLD / MDS-MLD)** | $< 5\%$ | $< 2\%$ | Significant dysplasia ($\ge 10\%$) in 1 (SLD) or $\ge 2$ (MLD) lineages |
| **MDS with Ring Sideroblasts (MDS-RS)** | $< 5\%$ | $< 2\%$ | $\ge 15\%$ ring sideroblasts, OR $\ge 5\%$ with *SF3B1* somatic mutation |
| **MDS with Increased Blasts 1 (MDS-IB1 / MDS-EB-1)** | $5.0\% - 9.9\%$ | $2.0\% - 4.9\%$ | Absence of Auer rods or AML-defining genetic lesions |
| **MDS with Increased Blasts 2 (MDS-IB2 / MDS-EB-2)** | $10.0\% - 19.9\%$ | $5.0\% - 19.9\%$ | Or presence of Auer rods regardless of blast count |
| **Acute Myeloid Leukemia (AML)** | $\ge 20\%$ | $\ge 20\%$ | Morphologic blast threshold for AML-NOS or AML-MR |
| **AML with Defining Recurrent Genetics** | Any / $\ge 10\%$ | Any / $\ge 10\%$ | *PML::RARA*, *RUNX1::RUNX1T1*, *CBFB::MYH11*, *KMT2A*, *MECOM*, *NPM1* |
| **Plasma Cell Myeloma** | $\ge 10\%$ | N/A | Bone marrow clonal plasma cells $\ge 10\%$ ($\ge 60\%$ is myeloma-defining biomarker) |

### 6. Reference Differential Ranges (Adult Bone Marrow Aspirate)

| Cell Type | Reference Percentage Range (%) |
|:---|:---:|
| **Myeloblasts** | 0.5 – 3.0% |
| **Promyelocytes** | 1.0 – 5.0% |
| **Myelocytes** | 5.0 – 15.0% |
| **Metamyelocytes** | 10.0 – 20.0% |
| **Band Neutrophils** | 10.0 – 20.0% |
| **Segmented Neutrophils** | 10.0 – 30.0% |
| **Eosinophils & Precursors** | 1.0 – 5.0% |
| **Basophils** | 0.1 – 1.0% |
| **Monocytes** | 1.0 – 4.0% |
| **Nucleated Erythroid Precursors (total)** | 15.0 – 30.0% |
| **Lymphocytes** | 5.0 – 15.0% |
| **Plasma Cells** | 0.5 – 3.0% |

---

## 💻 CLI Quickstart & Usage

The command-line interface provides single-case analysis, batch CSV processing, and built-in benchmark demos.

### 1. Batch CSV Processing
Process an entire cohort of aspirate differentials from CSV and write structured diagnostic sign-out reports:

```bash
# Process batch CSV and write to output file
python cli.py batch -i sample.csv -o out_results.csv

# Output directly to stdout
python cli.py batch -i sample.csv
```

### 2. Benchmark Demo Scenarios
Run pre-configured, validated clinical benchmarks:

```bash
python cli.py --demo normal      # Normocellular unremarkable marrow
python cli.py --demo aml         # Acute Myeloid Leukemia with NPM1/FLT3 mutations
python cli.py --demo mds_rs      # MDS with Ring Sideroblasts and SF3B1 mutation
python cli.py --demo aplastic    # Severe aplastic anemia pattern
python cli.py --demo all         # Run all benchmark scenarios
```

### 3. Direct Case Evaluation
Evaluate a case with specific cell counts and export as structured JSON or clinical text:

```bash
python cli.py \
  --case-id "BM-2026-0042" \
  --age 64 \
  --cellularity 85.0 \
  --blasts 125 \
  --promyelocytes 20 \
  --myelocytes 30 \
  --metamyelocytes 35 \
  --bands 40 \
  --segs 70 \
  --poly-normo 60 \
  --ortho-normo 40 \
  --lymphocytes 50 \
  --plasma-cells 10 \
  --json
```

### 4. Interactive Mode
Launch the step-by-step interactive prompt for aspirate counts:

```bash
python cli.py --interactive
```

---

## 🐍 Python API Quickstart

```python
from bone_marrow_differential import (
    BoneMarrowCellCounts,
    BoneMarrowDifferentialAnalyzer,
    ClinicalCaseInput,
    DysplasiaFeatures,
    format_clinical_report,
)

# 1. Enter 500-cell aspirate differential counts
counts = BoneMarrowCellCounts(
    blasts=35,              # 7.0% blasts
    promyelocytes=12,
    myelocytes=35,
    metamyelocytes=45,
    band_neutrophils=55,
    segmented_neutrophils=90,
    eosinophils=10,
    basophils=3,
    monocytes=10,
    pronormoblasts=8,
    basophilic_normoblasts=20,
    polychromatophilic_normoblasts=80,
    orthochromatophilic_normoblasts=45,
    lymphocytes=45,
    plasma_cells=7,
)

# 2. Build clinical case input
case = ClinicalCaseInput(
    case_id="BM-MDS-001",
    patient_age=68,
    counts=counts,
    core_cellularity_pct=65.0,
    peripheral_blood_blast_pct=1.5,
    dysplasia=DysplasiaFeatures(erythroid_dysplasia_pct=15.0),
)

# 3. Analyze against WHO 2022 / ICC criteria
report = BoneMarrowDifferentialAnalyzer.analyze(case)

# 4. Inspect report properties
print(f"Total Count: {report.total_cells_counted}")
print(f"Marrow Blasts: {report.marrow_blast_pct}%")
print(f"M:E Ratio: {report.me_ratio}:1")
print(f"Diagnosis: {report.primary_diagnostic_category}")
print(f"Subclass: {report.subclassification}")
print(f"IPSS-R Stratum: {report.ipss_r_blast_score_category}")

# 5. Format formatted sign-out report
print(format_clinical_report(report))
```

---

## 🧪 Testing & Verification

Run the comprehensive unit test suite:

```bash
python -m pytest -p no:zarr -v
```

Execute the batch CLI smoke verification:

```bash
python cli.py batch -i sample.csv -o out_smoke.csv
# Verify and clean up
python -c "import os; assert os.path.exists('out_smoke.csv'); os.remove('out_smoke.csv')"
```

---

## 📄 License

MIT License. See [LICENSE](LICENSE) for details.
