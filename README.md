# Bone Marrow Differential & Hematopathology Diagnostic Agent

An enterprise-grade, deterministic clinical decision-support and computational hematopathology engine implementing the **WHO 5th Edition (2022)** and **International Consensus Classification (ICC 2022)** criteria for bone marrow aspirate differentials, biopsy cellularity, dysplasia quantification, and myeloid/lymphoid neoplasm classification.

---

## Clinical & Scientific Overview

Manual 500-cell bone marrow aspirate evaluation remains the cornerstone of diagnostic hematopathology. This engine automates the validation, mathematical derivation, and guideline-driven classification of bone marrow findings.

### Key Clinical Capabilities

1. **500-Cell Aspirate Differential Analysis**:
   - Complete granular count across Granulocytic, Erythroid, Monocytic, Lymphoid, and Plasma cell lineages.
   - Total nucleated cell count validation and lineage percentage breakdown.
2. **Blast Percentage & Non-Erythroid Blast Percentage**:
   - $\text{Blast } \% = \frac{\text{Total Blasts}}{\text{Total Nucleated Cells}} \times 100$
   - $\text{Non-Erythroid Blast } \% = \frac{\text{Total Blasts}}{\text{Total Nucleated Cells} - \text{Total Erythroid Precursors}} \times 100$ (FAB criteria correlation).
3. **Myeloid-to-Erythroid (M:E) Ratio**:
   - Standard reference range: $1.5:1$ to $3.5:1$.
   - Identifies erythroid hyperplasia (hemolysis, blood loss, erythropoietin response) vs. granulocytic hyperplasia (infection, leukemoid reaction, CML).
4. **Age-Adjusted Core Biopsy Cellularity**:
   - Empirical standard formula: $\text{Expected Cellularity } (\%) = 100 - \text{Age}$.
   - Age-bounded reference ranges: Normocellular, Hypocellular, Severely Hypocellular ($<10\%$, Aplastic Anemia screening), and Hypercellular.
5. **Morphologic Dysplasia Quantification**:
   - Granular assessment across Erythroid, Granulocytic, and Megakaryocytic lineages.
   - Minimum threshold of $\ge 10\%$ dysplastic cells for single-lineage (SLD) or multilineage (MLD) dysplasia designation.
   - Ring sideroblast scoring with Perls' Prussian blue stain ($\ge 15\%$, or $\ge 5\%$ in the presence of somatic *SF3B1* mutations).
6. **WHO 2022 / ICC Diagnostic Classifier**:
   - Acute Myeloid Leukemia (AML) by blast count ($\ge 20\%$) and recurrent defining genetic lesions (*PML::RARA*, *RUNX1::RUNX1T1*, *CBFB::MYH11*, *KMT2A*, etc.).
   - Myelodysplastic Neoplasms (MDS-LB, MDS-IB1, MDS-IB2, MDS with isolated *del(5q)*, MDS-SF3B1).
   - Chronic Myelomonocytic Leukemia (CMML) and Plasma Cell Neoplasms / Multiple Myeloma ($\ge 10\%$ and $\ge 60\%$ biomarkers).
   - IPSS-R blast stratum scoring ($\le 2\%$, $>2-<5\%$, $5-10\%$, $>10\%$).

---

## Installation & Requirements

Requires **Python 3.9+** with standard library dependencies.

```bash
git clone https://github.com/abusuraihsakhri/bone-marrow-differential-agent.git
cd bone-marrow-differential-agent
```

---

## Command-Line Interface (CLI)

### 1. Run Pre-Configured Benchmark Scenarios

```bash
python cli.py --demo normal
python cli.py --demo aml
python cli.py --demo mds_rs
python cli.py --demo aplastic
```

### 2. Direct Argument Input with JSON Output

```bash
python cli.py --case-id BM-2026-089 --age 64 --cellularity 85 \
  --blasts 120 --segs 380 --json
```

### 3. Interactive Case Entry Mode

```bash
python cli.py --interactive
```

---

## Python API Usage

```python
from bone_marrow_differential import (
    BoneMarrowCellCounts,
    ClinicalCaseInput,
    DysplasiaFeatures,
    BoneMarrowDifferentialAnalyzer,
    format_clinical_report,
)

# Define 500-cell aspirate counts
counts = BoneMarrowCellCounts(
    blasts=15,
    promyelocytes=20,
    myelocytes=50,
    metamyelocytes=60,
    band_neutrophils=70,
    segmented_neutrophils=110,
    eosinophils=15,
    basophils=5,
    monocytes=15,
    pronormoblasts=5,
    basophilic_normoblasts=15,
    polychromatophilic_normoblasts=65,
    orthochromatophilic_normoblasts=35,
    lymphocytes=30,
    plasma_cells=5,
)

case = ClinicalCaseInput(
    case_id="PAT-9481",
    patient_age=68,
    counts=counts,
    core_cellularity_pct=65.0,
    peripheral_blood_blast_pct=1.0,
    dysplasia=DysplasiaFeatures(
        erythroid_dysplasia_pct=15.0,
        ring_sideroblasts_pct=18.0,
        sf3b1_mutation_detected=True
    )
)

report = BoneMarrowDifferentialAnalyzer.analyze(case)
print(format_clinical_report(report))
```

---

## Verification & Unit Test Suite

Run the full unit test suite with 31 automated test cases:

```bash
python -m unittest test_bone_marrow_differential.py -v
```

---

## License

MIT License. Authored and maintained by Dr. Abu Suraih Sakhri.
