# Bone Marrow Differential Agent

> **Domain:** Digital Pathology & Quantitative Histopathology  
> **Reference Guidelines & Standards:** `College of American Pathologists (CAP) Synoptic Protocols & DICOM WSI`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**Bone Marrow Differential Agent** is an advanced analytical and computational platform implementing 500-Cell Aspirate Differential & Blast Fraction Harmonizer.

Bone Marrow Differential & Hematopathology Analysis Engine
==========================================================
Comprehensive hematopathology system for 500-cell aspirate differential counts,
Myeloid-to-Erythroid (M:E) ratio computation, age-adjusted core biopsy cellularity,
Perls' Prussian blue iron stores/ring sideroblast quantification, dysplasia scoring,
and WHO 2022 / ICC hematologic neoplasm diagnostic classification.

Standards & Criteria:
- WHO Classification of Haematolymphoid Tumours (5th Edition, 2022)
- International Consensus Classification (ICC 2022)
- International Prognostic Scoring System (IPSS-R / IPSS-M) Blast Strata
- Standard Hematopathology Manual Differential Count Benchmarks (500-cell standard)

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`BlastMeasurements`** — dedicated module for blast measurements evaluation and state verification.
- **`ConsistencyVerdict`** — dedicated module for consistency verdict evaluation and state verification.
- **`Lineage`** — dedicated module for lineage evaluation and state verification.
- **`CellularityStatus`** — dedicated module for cellularity status evaluation and state verification.
- **`DysplasiaDegree`** — dedicated module for dysplasia degree evaluation and state verification.
- **`IronStoreGrade`** — dedicated module for iron store grade evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  return (
  Formula: Blasts / (Total Cells - Total Erythroid) * 100
  """Calculates percentage for every counted cell category."""
  """Calculates Myeloid-to-Erythroid ratio."""
  weighted_score = weight * (pct / 100.0)
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --interactive <value> --demo <value> --case-id <value> --age <value>
```

### Parameter Reference
- `--interactive`: Specifies input measurement or parameter value.
- `--demo`: Specifies input measurement or parameter value.
- `--case-id`: Specifies input measurement or parameter value.
- `--age`: Specifies input measurement or parameter value.
- `--cellularity`: Specifies input measurement or parameter value.
- `--pb-blasts`: Specifies input measurement or parameter value.
- `--pb-monos`: Specifies input measurement or parameter value.
- `--blasts`: Specifies input measurement or parameter value.
- `--promyelocytes`: Specifies input measurement or parameter value.
- `--myelocytes`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `case_id` | Parameter / observation metric | Required |
| `patient_synthetic_id` | Parameter / observation metric | Required |
| `metric_primary` | Parameter / observation metric | Required |
| `metric_secondary` | Parameter / observation metric | Required |
| `is_stat` | Parameter / observation metric | Required |
| `status_flag` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t bone-marrow-differential-agent .
docker run -p 8000:8000 bone-marrow-differential-agent
```
