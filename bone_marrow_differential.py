"""
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
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Dict, Any, List, Optional, Tuple, Union


class Lineage(str, Enum):
    GRANULOCYTIC = "Granulocytic / Myeloid"
    ERYTHROID = "Erythroid"
    MONOCYTIC = "Monocytic"
    LYMPHOID = "Lymphoid"
    MEGAKARYOCYTIC = "Megakaryocytic"
    PLASMA_CELL = "Plasma Cell"
    OTHER = "Other / Histiocytic / Mast"


class CellularityStatus(str, Enum):
    SEVERELY_HYPOCELLULAR = "Severely Hypocellular"
    HYPOCELLULAR = "Hypocellular"
    NORMOCELLULAR = "Normocellular"
    HYPERCELLULAR = "Hypercellular"
    SEVERELY_HYPERCELLULAR = "Severely Hypercellular"


class DysplasiaDegree(str, Enum):
    NONE = "No Significant Dysplasia (<10%)"
    SINGLE_LINEAGE = "Single Lineage Dysplasia (SLD, >=10% in 1 lineage)"
    MULTILINEAGE = "Multilineage Dysplasia (MLD, >=10% in >=2 lineages)"


class IronStoreGrade(int, Enum):
    GRADE_0 = 0  # None / Absent
    GRADE_1 = 1  # Severely decreased
    GRADE_2 = 2  # Slightly decreased
    GRADE_3 = 3  # Normal
    GRADE_4 = 4  # Moderately increased
    GRADE_5 = 5  # Markedly increased
    GRADE_6 = 6  # Very large clumps / massive iron overload


@dataclass
class BoneMarrowCellCounts:
    """
    Detailed manual differential count of nucleated bone marrow cells.
    Recommended standard count is 500 nucleated cells from aspirate smears.
    """
    # Granulocytic / Myeloid lineage
    blasts: int = 0
    promyelocytes: int = 0
    myelocytes: int = 0
    metamyelocytes: int = 0
    band_neutrophils: int = 0
    segmented_neutrophils: int = 0
    eosinophils: int = 0
    basophils: int = 0

    # Monocytic lineage
    monocytes: int = 0

    # Erythroid lineage
    pronormoblasts: int = 0
    basophilic_normoblasts: int = 0
    polychromatophilic_normoblasts: int = 0
    orthochromatophilic_normoblasts: int = 0

    # Lymphoid & Plasma cells
    lymphocytes: int = 0
    plasma_cells: int = 0

    # Other lineages
    megakaryocytes: int = 0
    histiocytes: int = 0
    mast_cells: int = 0
    other_cells: int = 0

    def total_count(self) -> int:
        """Returns total nucleated cell count."""
        return (
            self.blasts + self.promyelocytes + self.myelocytes + self.metamyelocytes +
            self.band_neutrophils + self.segmented_neutrophils + self.eosinophils +
            self.basophils + self.monocytes + self.pronormoblasts +
            self.basophilic_normoblasts + self.polychromatophilic_normoblasts +
            self.orthochromatophilic_normoblasts + self.lymphocytes +
            self.plasma_cells + self.megakaryocytes + self.histiocytes +
            self.mast_cells + self.other_cells
        )

    def total_granulocytic(self) -> int:
        """
        Total myeloid/granulocytic cells including precursors:
        Blasts (myeloid), Promyelocytes, Myelocytes, Metamyelocytes, Bands, Segs, Eosinophils, Basophils.
        """
        return (
            self.blasts + self.promyelocytes + self.myelocytes + self.metamyelocytes +
            self.band_neutrophils + self.segmented_neutrophils + self.eosinophils +
            self.basophils
        )

    def total_erythroid(self) -> int:
        """Total nucleated erythroid precursors."""
        return (
            self.pronormoblasts + self.basophilic_normoblasts +
            self.polychromatophilic_normoblasts + self.orthochromatophilic_normoblasts
        )

    def total_monocytic(self) -> int:
        return self.monocytes

    def total_lymphoid(self) -> int:
        return self.lymphocytes + self.plasma_cells

    def blast_percentage(self) -> float:
        """Blast percentage of total nucleated bone marrow cells."""
        tot = self.total_count()
        if tot == 0:
            return 0.0
        return round((self.blasts / tot) * 100.0, 2)

    def non_erythroid_blast_percentage(self) -> float:
        """
        Blast percentage of non-erythroid nucleated cells (FAB criteria).
        Formula: Blasts / (Total Cells - Total Erythroid) * 100
        """
        non_erythroid_total = self.total_count() - self.total_erythroid()
        if non_erythroid_total <= 0:
            return 0.0
        return round((self.blasts / non_erythroid_total) * 100.0, 2)

    def myeloid_to_erythroid_ratio(self) -> float:
        """
        Myeloid to Erythroid (M:E) Ratio.
        Normal adult reference range is typically 1.5:1 to 3.5:1 (or 2.0 to 4.0:1).
        """
        erythroid = self.total_erythroid()
        if erythroid == 0:
            return float("inf") if self.total_granulocytic() > 0 else 0.0
        return round(self.total_granulocytic() / erythroid, 2)

    def percentages(self) -> Dict[str, float]:
        """Calculates percentage for every counted cell category."""
        tot = self.total_count()
        if tot == 0:
            return {k: 0.0 for k in self.__dict__.keys()}
        return {
            k: round((v / tot) * 100.0, 2)
            for k, v in self.__dict__.items()
            if isinstance(v, (int, float))
        }

    def validate(self) -> List[str]:
        """Validates count entries for clinical consistency and adequacy."""
        issues = []
        for k, v in self.__dict__.items():
            if isinstance(v, (int, float)) and v < 0:
                issues.append(f"Negative count not permitted for {k}: {v}")
        tot = self.total_count()
        if tot == 0:
            issues.append("Total nucleated cell count is 0; at least 100 cells required for differential.")
        elif tot < 200:
            issues.append(f"Count of {tot} cells is below standard minimal threshold (200-500 cells recommended).")
        elif tot < 500:
            issues.append(f"Count of {tot} cells is acceptable but below optimal 500-cell aspirate standard.")
        return issues


@dataclass
class CellularityAssessment:
    """Age-adjusted core biopsy cellularity analysis."""
    patient_age: int
    observed_cellularity_pct: float
    expected_cellularity_pct: float = field(init=False)
    lower_normal_limit_pct: float = field(init=False)
    upper_normal_limit_pct: float = field(init=False)
    status: CellularityStatus = field(init=False)
    interpretation: str = field(init=False)

    def __post_init__(self):
        if not (0 <= self.patient_age <= 125):
            raise ValueError(f"Patient age must be between 0 and 125, got {self.patient_age}")
        if not (0.0 <= self.observed_cellularity_pct <= 100.0):
            raise ValueError(f"Cellularity percentage must be between 0 and 100, got {self.observed_cellularity_pct}")

        # Standard rule: Expected cellularity (%) = 100 - age
        # Pediatric adjustment: neonates ~ 100%, young infants ~ 80-90%
        if self.patient_age < 2:
            self.expected_cellularity_pct = 95.0
        elif self.patient_age < 10:
            self.expected_cellularity_pct = 85.0
        else:
            self.expected_cellularity_pct = max(10.0, float(100 - self.patient_age))

        # Normal standard window is expected +/- 15% (bounded between 5% and 95%)
        self.lower_normal_limit_pct = max(5.0, self.expected_cellularity_pct - 15.0)
        self.upper_normal_limit_pct = min(95.0, self.expected_cellularity_pct + 15.0)

        # Classification
        if self.observed_cellularity_pct < 10.0:
            self.status = CellularityStatus.SEVERELY_HYPOCELLULAR
            self.interpretation = f"Severe hypocellularity ({self.observed_cellularity_pct}%). Evaluate for aplastic anemia or profound hypoplasia."
        elif self.observed_cellularity_pct < self.lower_normal_limit_pct:
            self.status = CellularityStatus.HYPOCELLULAR
            self.interpretation = f"Hypocellular for age {self.patient_age} (observed: {self.observed_cellularity_pct}%, expected: {self.expected_cellularity_pct}%)."
        elif self.observed_cellularity_pct > min(95.0, self.upper_normal_limit_pct + 15.0):
            self.status = CellularityStatus.SEVERELY_HYPERCELLULAR
            self.interpretation = f"Marked hypercellularity ({self.observed_cellularity_pct}%). High suspicion for marrow proliferation or leukemic infiltration."
        elif self.observed_cellularity_pct > self.upper_normal_limit_pct:
            self.status = CellularityStatus.HYPERCELLULAR
            self.interpretation = f"Hypercellular for age {self.patient_age} (observed: {self.observed_cellularity_pct}%, expected: {self.expected_cellularity_pct}%)."
        else:
            self.status = CellularityStatus.NORMOCELLULAR
            self.interpretation = f"Normocellular marrow for age {self.patient_age} (observed {self.observed_cellularity_pct}% within {self.lower_normal_limit_pct:.0f}-{self.upper_normal_limit_pct:.0f}%)."


@dataclass
class DysplasiaFeatures:
    """Dysplasia assessment across 3 main lineages."""
    erythroid_dysplasia_pct: float = 0.0  # Nuclear budding, multinuclearity, ring sideroblasts, megaloblastoid
    granulocytic_dysplasia_pct: float = 0.0  # Hypogranulation, pseudo-Pelger-Huet, hypersegmentation
    megakaryocytic_dysplasia_pct: float = 0.0  # Micromegakaryocytes, multinucleated/separated lobes

    # Dysplastic qualitative findings
    auer_rods_present: bool = False
    ring_sideroblasts_pct: float = 0.0  # % of erythroid precursors
    sf3b1_mutation_detected: bool = False

    def dysplastic_lineages_count(self) -> int:
        """Significant dysplasia requires >= 10% morphologically abnormal cells in that lineage."""
        count = 0
        if self.erythroid_dysplasia_pct >= 10.0 or (self.ring_sideroblasts_pct >= 15.0 or (self.sf3b1_mutation_detected and self.ring_sideroblasts_pct >= 5.0)):
            count += 1
        if self.granulocytic_dysplasia_pct >= 10.0:
            count += 1
        if self.megakaryocytic_dysplasia_pct >= 10.0:
            count += 1
        return count

    def get_dysplasia_degree(self) -> DysplasiaDegree:
        cnt = self.dysplastic_lineages_count()
        if cnt >= 2:
            return DysplasiaDegree.MULTILINEAGE
        elif cnt == 1:
            return DysplasiaDegree.SINGLE_LINEAGE
        return DysplasiaDegree.NONE


@dataclass
class ClinicalCaseInput:
    """Full clinical case input parameters for bone marrow evaluation."""
    case_id: str
    patient_age: int
    counts: BoneMarrowCellCounts
    core_cellularity_pct: float = 50.0
    peripheral_blood_blast_pct: float = 0.0
    peripheral_blood_monocyte_abs_k_ul: float = 0.5  # x10^9/L or k/uL
    dysplasia: DysplasiaFeatures = field(default_factory=DysplasiaFeatures)
    iron_store_grade: Optional[IronStoreGrade] = IronStoreGrade.GRADE_3
    cytogenetics_or_mutations: List[str] = field(default_factory=list)
    flow_cytometry_markers: Dict[str, str] = field(default_factory=dict)
    clinical_history: str = ""


@dataclass
class BoneMarrowReport:
    """Comprehensive diagnostic and analytical report."""
    case_id: str
    patient_age: int
    total_cells_counted: int
    marrow_blast_pct: float
    non_erythroid_blast_pct: float
    me_ratio: float
    cellularity: CellularityAssessment
    dysplasia_degree: DysplasiaDegree
    dysplastic_lineages_count: int
    primary_diagnostic_category: str
    subclassification: str
    who_2022_criteria_matched: List[str]
    ipss_r_blast_score_category: str
    critical_alerts: List[str]
    advisory_recommendations: List[str]
    differential_percentages: Dict[str, float]

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["cellularity"]["status"] = self.cellularity.status.value
        d["dysplasia_degree"] = self.dysplasia_degree.value
        return d

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


class BoneMarrowDifferentialAnalyzer:
    """
    Expert Analyzer implementing WHO 5th Ed (2022), ICC 2022, and IPSS-R hematopathology guidelines.
    """

    # High-risk AML-defining cytogenetic abnormalities (WHO 2022 defines AML regardless of blast count >= 20%)
    AML_DEFINING_MUTATIONS = [
        "PML::RARA", "PML-RARA", "t(15;17)",
        "RUNX1::RUNX1T1", "RUNX1-RUNX1T1", "t(8;21)",
        "CBFB::MYH11", "CBFB-MYH11", "inv(16)", "t(16;16)",
        "KMT2A", "t(9;11)", "MECOM", "inv(3)", "t(3;3)",
        "NUP98", "NPM1"
    ]

    @staticmethod
    def calculate_me_ratio(granulocytic_count: int, erythroid_count: int) -> float:
        """Calculates Myeloid-to-Erythroid ratio."""
        if erythroid_count == 0:
            return float("inf") if granulocytic_count > 0 else 0.0
        return round(granulocytic_count / erythroid_count, 2)

    @classmethod
    def analyze(cls, case: ClinicalCaseInput) -> BoneMarrowReport:
        """Executes full diagnostic and numerical synthesis on bone marrow data."""
        validation_issues = case.counts.validate()
        total_counted = case.counts.total_count()
        blast_pct = case.counts.blast_percentage()
        non_erythroid_blast_pct = case.counts.non_erythroid_blast_percentage()
        me_ratio = case.counts.myeloid_to_erythroid_ratio()
        percentages = case.counts.percentages()

        cellularity = CellularityAssessment(
            patient_age=case.patient_age,
            observed_cellularity_pct=case.core_cellularity_pct
        )

        dysp = case.dysplasia
        dysplasia_degree = dysp.get_dysplasia_degree()
        dysp_count = dysp.dysplastic_lineages_count()

        # Check for AML-defining genetic abnormalities
        aml_defining_genetics = []
        for mut in case.cytogenetics_or_mutations:
            for aml_def in cls.AML_DEFINING_MUTATIONS:
                if aml_def.lower() in mut.lower():
                    aml_defining_genetics.append(mut)

        criteria_matched = []
        critical_alerts = []
        advisories = []

        # Check validation warnings
        for issue in validation_issues:
            advisories.append(f"Validation Notice: {issue}")

        # Determine Primary Category & Subclassification
        primary_diag = "Unclassified / Reactive"
        subclass = "Pending comprehensive correlation"

        # 1. ACUTE MYELOID LEUKEMIA (AML)
        if blast_pct >= 20.0 or case.peripheral_blood_blast_pct >= 20.0:
            primary_diag = "Acute Myeloid Leukemia (AML)"
            criteria_matched.append(f"Blasts >= 20% (Bone Marrow: {blast_pct}%, PB: {case.peripheral_blood_blast_pct}%)")
            critical_alerts.append("CRITICAL: Blasts >= 20% confirms Acute Leukemia / AML threshold.")
            if aml_defining_genetics:
                subclass = f"AML with recurrent genetic abnormalities ({', '.join(aml_defining_genetics)})"
                criteria_matched.append(f"Recurrent genetics: {aml_defining_genetics}")
            elif dysp_count >= 2:
                subclass = "AML with Myelodysplasia-Related Gene Mutations / Cytogenetics (AML-MR)"
                criteria_matched.append(f"Multilineage dysplasia ({dysp_count} lineages) with blast crisis")
            else:
                subclass = "AML, Not Otherwise Specified (AML-NOS)"

        elif aml_defining_genetics and (blast_pct >= 10.0 or case.peripheral_blood_blast_pct >= 10.0 or any("PML" in g or "t(15;17)" in g for g in aml_defining_genetics)):
            primary_diag = "Acute Myeloid Leukemia (AML) by Defining Genetics"
            subclass = f"AML defined by recurrent genetic abnormality: {', '.join(aml_defining_genetics)}"
            criteria_matched.append(f"WHO 2022 / ICC: AML diagnosed with recurrent abnormality ({aml_defining_genetics}) despite blasts {blast_pct}%")
            critical_alerts.append("CRITICAL: AML-defining genetic lesion identified.")

        # 2. MYELODYSPLASTIC NEOPLASMS (MDS) / MDS-IB / MDS-LB
        elif (5.0 <= blast_pct < 20.0) or (2.0 <= case.peripheral_blood_blast_pct < 20.0) or dysp.auer_rods_present:
            primary_diag = "Myelodysplastic Neoplasm (MDS) with Increased Blasts"
            if blast_pct >= 10.0 or case.peripheral_blood_blast_pct >= 5.0 or dysp.auer_rods_present:
                subclass = "MDS with Increased Blasts 2 (MDS-IB2 / ICC: MDS with Excess Blasts 2)"
                criteria_matched.append(f"Blasts 10-19% in marrow ({blast_pct}%) or 5-19% in PB ({case.peripheral_blood_blast_pct}%) or Auer rods present")
            else:
                subclass = "MDS with Increased Blasts 1 (MDS-IB1 / ICC: MDS with Excess Blasts 1)"
                criteria_matched.append(f"Blasts 5-9% in marrow ({blast_pct}%) or 2-4% in PB ({case.peripheral_blood_blast_pct}%)")
            critical_alerts.append(f"ELEVATED BLASTS: {blast_pct}% marrow blasts indicates high-risk MDS/advanced myeloid neoplasm.")

        elif dysp_count >= 1:
            primary_diag = "Myelodysplastic Neoplasm (MDS)"
            # Check for isolated del(5q)
            is_del5q = any("5q" in g.lower() or "del(5q)" in g.lower() for g in case.cytogenetics_or_mutations)
            is_sf3b1 = dysp.sf3b1_mutation_detected or any("sf3b1" in g.lower() for g in case.cytogenetics_or_mutations)

            if is_del5q and blast_pct < 5.0 and case.peripheral_blood_blast_pct < 2.0:
                subclass = "MDS with isolated del(5q)"
                criteria_matched.append("del(5q) cytogenetics with blasts < 5%")
            elif (dysp.ring_sideroblasts_pct >= 15.0 or (is_sf3b1 and dysp.ring_sideroblasts_pct >= 5.0)) and blast_pct < 5.0:
                subclass = "MDS with Ring Sideroblasts (MDS-RS / MDS-SF3B1)"
                criteria_matched.append(f"Ring sideroblasts {dysp.ring_sideroblasts_pct}% (SF3B1: {is_sf3b1}) with blasts < 5%")
            elif dysp_count >= 2:
                subclass = "MDS with Multilineage Dysplasia (MDS-MLD / MDS-LB with MLD)"
                criteria_matched.append(f"Dysplasia >=10% in {dysp_count} lineages (Erythroid: {dysp.erythroid_dysplasia_pct}%, Granulocytic: {dysp.granulocytic_dysplasia_pct}%, Megakaryocytic: {dysp.megakaryocytic_dysplasia_pct}%)")
            else:
                subclass = "MDS with Single Lineage Dysplasia (MDS-SLD / MDS-LB with SLD)"
                criteria_matched.append("Dysplasia >=10% restricted to 1 lineage with blasts < 5%")

        # 3. CHRONIC MYELOMONOCYTIC LEUKEMIA (CMML)
        elif (case.counts.monocytes / max(1, total_counted) >= 0.10 or case.peripheral_blood_monocyte_abs_k_ul >= 1.0) and blast_pct < 20.0:
            if case.counts.monocytes / max(1, total_counted) >= 0.10:
                primary_diag = "Myelodysplastic / Myeloproliferative Neoplasm (MDS/MPN)"
                subclass = "Chronic Myelomonocytic Leukemia (CMML) Pattern"
                criteria_matched.append(f"Monocytosis >= 10% in marrow differential ({percentages.get('monocytes', 0)}%)")
                advisories.append("Correlate with absolute persistent peripheral blood monocytosis >= 1.0 x 10^9/L.")

        # 4. PLASMA CELL NEOPLASMS
        elif percentages.get("plasma_cells", 0.0) >= 10.0:
            primary_diag = "Plasma Cell Neoplasm"
            pct_plasma = percentages.get("plasma_cells", 0.0)
            if pct_plasma >= 60.0:
                subclass = "Multiple Myeloma (Myeloma-defining biomarker: Plasma cells >= 60%)"
                criteria_matched.append(f"Marrow plasmacytosis >= 60% ({pct_plasma}%)")
                critical_alerts.append("CRITICAL: Plasma cells >= 60% meets SLiM-CRAB criteria for Multiple Myeloma.")
            else:
                subclass = f"Plasma Cell Myeloma / MGUS (Plasma cells: {pct_plasma}%)"
                criteria_matched.append(f"Bone marrow plasma cells {pct_plasma}% (threshold >= 10%)")
                advisories.append("Evaluate CRAB criteria (Hypercalcemia, Renal failure, Anemia, Bone lytic lesions) and Serum Free Light Chains.")

        # 5. REACTIVE / PHYSIOLOGICAL / APLASTIC PATTERNS
        elif cellularity.status == CellularityStatus.SEVERELY_HYPOCELLULAR and blast_pct < 5.0 and dysp_count == 0:
            primary_diag = "Bone Marrow Failure Syndrome"
            subclass = "Aplastic Anemia / Severe Hypoplasia"
            criteria_matched.append(f"Severe hypocellularity ({cellularity.observed_cellularity_pct}%) without increased blasts or dysplasia")
            critical_alerts.append("CRITICAL: Severe marrow hypocellularity. Exclude PNH clone and Fanconi anemia.")

        elif me_ratio > 4.5:
            primary_diag = "Myeloid Hyperplasia / Shift"
            subclass = f"Granulocytic Hyperplasia (M:E Ratio: {me_ratio}:1, normal 1.5-3.5:1)"
            criteria_matched.append(f"Elevated M:E ratio ({me_ratio}:1) with normal blast percentage ({blast_pct}%)")
            advisories.append("Assess for underlying leukemoid reaction, systemic bacterial infection, G-CSF administration, or early CML (test BCR-ABL1).")

        elif me_ratio < 1.2 and me_ratio > 0:
            primary_diag = "Erythroid Hyperplasia"
            subclass = f"Erythroid Hyperplasia (M:E Ratio: {me_ratio}:1, inverted)"
            criteria_matched.append(f"Inverted/decreased M:E ratio ({me_ratio}:1) showing erythroid predominance")
            advisories.append("Assess for hemolytic anemia, blood loss recovery, erythropoietin therapy, or thalassemia.")

        else:
            primary_diag = "Normocellular / Unremarkable Bone Marrow"
            subclass = "Morphologically within normal age-adjusted reference ranges"
            criteria_matched.append(f"Normocellular marrow, normal M:E ratio ({me_ratio}:1), blasts < 5% ({blast_pct}%)")

        # IPSS-R Blast Category
        if blast_pct <= 2.0:
            ipss_r_blast = "<= 2% (IPSS-R Score: 0)"
        elif blast_pct > 2.0 and blast_pct < 5.0:
            ipss_r_blast = "> 2% to < 5% (IPSS-R Score: 1)"
        elif 5.0 <= blast_pct <= 10.0:
            ipss_r_blast = "5% to 10% (IPSS-R Score: 2)"
        else:
            ipss_r_blast = "> 10% (IPSS-R Score: 3)"

        # Iron stain correlation
        if case.iron_store_grade is not None:
            if case.iron_store_grade in (IronStoreGrade.GRADE_0, IronStoreGrade.GRADE_1):
                advisories.append(f"Iron Stores: Grade {case.iron_store_grade.value} (Depleted/Severely Decreased). Suggests iron deficiency.")
            elif case.iron_store_grade in (IronStoreGrade.GRADE_5, IronStoreGrade.GRADE_6):
                advisories.append(f"Iron Stores: Grade {case.iron_store_grade.value} (Markedly Increased). Correlate with transfusion history / hemochromatosis.")

        return BoneMarrowReport(
            case_id=case.case_id,
            patient_age=case.patient_age,
            total_cells_counted=total_counted,
            marrow_blast_pct=blast_pct,
            non_erythroid_blast_pct=non_erythroid_blast_pct,
            me_ratio=me_ratio,
            cellularity=cellularity,
            dysplasia_degree=dysplasia_degree,
            dysplastic_lineages_count=dysp_count,
            primary_diagnostic_category=primary_diag,
            subclassification=subclass,
            who_2022_criteria_matched=criteria_matched,
            ipss_r_blast_score_category=ipss_r_blast,
            critical_alerts=critical_alerts,
            advisory_recommendations=advisories,
            differential_percentages=percentages
        )


def format_clinical_report(report: BoneMarrowReport) -> str:
    """Formats the BoneMarrowReport into a clean, legible clinical text summary."""
    lines = []
    lines.append("=" * 78)
    lines.append(f" BONE MARROW DIFFERENTIAL & HISTOPATHOLOGY REPORT : {report.case_id}")
    lines.append("=" * 78)
    lines.append(f"Patient Age: {report.patient_age} yrs | Total Counted: {report.total_cells_counted} cells")
    lines.append(f"Marrow Blast %: {report.marrow_blast_pct:.2f}% | Non-Erythroid Blast %: {report.non_erythroid_blast_pct:.2f}%")
    lines.append(f"Myeloid:Erythroid (M:E) Ratio: {report.me_ratio:.2f}:1")
    lines.append(f"Core Cellularity: {report.cellularity.observed_cellularity_pct:.1f}% (Expected: {report.cellularity.expected_cellularity_pct:.1f}%) -> {report.cellularity.status.value}")
    lines.append(f"Dysplasia Assessment: {report.dysplasia_degree.value} ({report.dysplastic_lineages_count} lineages)")
    lines.append(f"IPSS-R Blast Stratum: {report.ipss_r_blast_score_category}")
    lines.append("-" * 78)
    lines.append(f"PRIMARY DIAGNOSIS: {report.primary_diagnostic_category}")
    lines.append(f"SUBCLASSIFICATION: {report.subclassification}")
    lines.append("-" * 78)

    if report.who_2022_criteria_matched:
        lines.append("Diagnostic Criteria Matched:")
        for crit in report.who_2022_criteria_matched:
            lines.append(f"  - {crit}")

    if report.critical_alerts:
        lines.append("\n[!] CRITICAL PATHOLOGY ALERTS:")
        for alert in report.critical_alerts:
            lines.append(f"  * {alert}")

    if report.advisory_recommendations:
        lines.append("\n[*] CLINICAL ADVISORIES & RECOMMENDATIONS:")
        for adv in report.advisory_recommendations:
            lines.append(f"  * {adv}")

    lines.append("\nCELL DIFFERENTIAL BREAKDOWN (%):")
    for cell_type, pct in report.differential_percentages.items():
        if pct > 0:
            lines.append(f"  - {cell_type.replace('_', ' ').title():<32}: {pct:>6.2f}%")

    lines.append("=" * 78)
    return "\n".join(lines)
