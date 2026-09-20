"""Bone marrow differential calculations with conservative interpretation flags.

The module is an interpretation aid, not a standalone diagnostic or treatment
system. Disease classification requires integrated clinical, blood, morphology,
flow, cytogenetic and/or molecular data that are not fully represented here.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


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
    SINGLE_LINEAGE = "Single Lineage Dysplasia (>=10% in 1 lineage)"
    MULTILINEAGE = "Multilineage Dysplasia (>=10% in >=2 lineages)"


class IronStoreGrade(int, Enum):
    GRADE_0 = 0
    GRADE_1 = 1
    GRADE_2 = 2
    GRADE_3 = 3
    GRADE_4 = 4
    GRADE_5 = 5
    GRADE_6 = 6


@dataclass
class BoneMarrowCellCounts:
    blasts: int = 0
    promyelocytes: int = 0
    myelocytes: int = 0
    metamyelocytes: int = 0
    band_neutrophils: int = 0
    segmented_neutrophils: int = 0
    eosinophils: int = 0
    basophils: int = 0
    monocytes: int = 0
    pronormoblasts: int = 0
    basophilic_normoblasts: int = 0
    polychromatophilic_normoblasts: int = 0
    orthochromatophilic_normoblasts: int = 0
    lymphocytes: int = 0
    plasma_cells: int = 0
    megakaryocytes: int = 0
    histiocytes: int = 0
    mast_cells: int = 0
    other_cells: int = 0

    def total_count(self) -> int:
        return sum(self.__dict__.values())

    def total_granulocytic(self) -> int:
        return sum(
            getattr(self, name)
            for name in (
                "blasts", "promyelocytes", "myelocytes", "metamyelocytes",
                "band_neutrophils", "segmented_neutrophils", "eosinophils", "basophils",
            )
        )

    def total_erythroid(self) -> int:
        return sum(
            getattr(self, name)
            for name in (
                "pronormoblasts", "basophilic_normoblasts",
                "polychromatophilic_normoblasts", "orthochromatophilic_normoblasts",
            )
        )

    def total_monocytic(self) -> int:
        return self.monocytes

    def total_lymphoid(self) -> int:
        return self.lymphocytes + self.plasma_cells

    def blast_percentage(self) -> float:
        total = self.total_count()
        return round(self.blasts / total * 100, 2) if total else 0.0

    def non_erythroid_blast_percentage(self) -> float:
        denominator = self.total_count() - self.total_erythroid()
        return round(self.blasts / denominator * 100, 2) if denominator > 0 else 0.0

    def myeloid_to_erythroid_ratio(self) -> float:
        erythroid = self.total_erythroid()
        if erythroid == 0:
            return float("inf") if self.total_granulocytic() else 0.0
        return round(self.total_granulocytic() / erythroid, 2)

    def percentages(self) -> Dict[str, float]:
        total = self.total_count()
        if not total:
            return {name: 0.0 for name in self.__dict__}
        return {name: round(value / total * 100, 2) for name, value in self.__dict__.items()}

    def validate(self) -> List[str]:
        issues: List[str] = []
        for name, value in self.__dict__.items():
            if isinstance(value, bool) or not isinstance(value, int):
                issues.append(f"Count for {name} must be an integer, got {value!r}")
            elif value < 0:
                issues.append(f"Negative count not permitted for {name}: {value}")
        total = self.total_count()
        if total == 0:
            issues.append("Total nucleated cell count is 0; enter a marrow differential before analysis.")
        elif total < 200:
            issues.append(
                f"Count of {total} cells is below the usual minimum for a reliable marrow differential; "
                "500 nucleated cells are preferred when material permits."
            )
        elif total < 500:
            issues.append(f"Count of {total} cells is below the preferred 500-cell marrow differential standard.")
        return issues


@dataclass
class CellularityAssessment:
    patient_age: int
    observed_cellularity_pct: float
    expected_cellularity_pct: float = field(init=False)
    lower_normal_limit_pct: float = field(init=False)
    upper_normal_limit_pct: float = field(init=False)
    status: CellularityStatus = field(init=False)
    interpretation: str = field(init=False)

    def __post_init__(self) -> None:
        if not 0 <= self.patient_age <= 125:
            raise ValueError(f"Patient age must be between 0 and 125, got {self.patient_age}")
        if not 0 <= self.observed_cellularity_pct <= 100:
            raise ValueError("Cellularity percentage must be between 0 and 100")
        if self.patient_age < 2:
            expected = 95.0
        elif self.patient_age < 10:
            expected = 85.0
        else:
            expected = max(10.0, float(100 - self.patient_age))
        self.expected_cellularity_pct = expected
        self.lower_normal_limit_pct = max(5.0, expected - 15)
        self.upper_normal_limit_pct = min(95.0, expected + 15)
        observed = self.observed_cellularity_pct
        if observed < 10:
            self.status = CellularityStatus.SEVERELY_HYPOCELLULAR
        elif observed < self.lower_normal_limit_pct:
            self.status = CellularityStatus.HYPOCELLULAR
        elif observed > min(95, self.upper_normal_limit_pct + 15):
            self.status = CellularityStatus.SEVERELY_HYPERCELLULAR
        elif observed > self.upper_normal_limit_pct:
            self.status = CellularityStatus.HYPERCELLULAR
        else:
            self.status = CellularityStatus.NORMOCELLULAR
        self.interpretation = (
            f"Observed cellularity {observed:.1f}% is {self.status.value.lower()} relative to an "
            f"approximate age-based reference of {expected:.1f}% "
            f"({self.lower_normal_limit_pct:.0f}-{self.upper_normal_limit_pct:.0f}% reference window)."
        )


@dataclass
class DysplasiaFeatures:
    erythroid_dysplasia_pct: float = 0.0
    granulocytic_dysplasia_pct: float = 0.0
    megakaryocytic_dysplasia_pct: float = 0.0
    auer_rods_present: bool = False
    ring_sideroblasts_pct: float = 0.0
    sf3b1_mutation_detected: bool = False

    def validate(self) -> List[str]:
        issues = []
        for name in (
            "erythroid_dysplasia_pct", "granulocytic_dysplasia_pct",
            "megakaryocytic_dysplasia_pct", "ring_sideroblasts_pct",
        ):
            value = getattr(self, name)
            if not 0 <= value <= 100:
                issues.append(f"{name} must be between 0 and 100, got {value}")
        return issues

    def dysplastic_lineages_count(self) -> int:
        return sum(
            value >= 10
            for value in (
                self.erythroid_dysplasia_pct,
                self.granulocytic_dysplasia_pct,
                self.megakaryocytic_dysplasia_pct,
            )
        )

    def get_dysplasia_degree(self) -> DysplasiaDegree:
        count = self.dysplastic_lineages_count()
        return (
            DysplasiaDegree.MULTILINEAGE if count >= 2
            else DysplasiaDegree.SINGLE_LINEAGE if count == 1
            else DysplasiaDegree.NONE
        )


@dataclass
class ClinicalCaseInput:
    case_id: str
    patient_age: int
    counts: BoneMarrowCellCounts
    core_cellularity_pct: float = 50.0
    peripheral_blood_blast_pct: float = 0.0
    peripheral_blood_monocyte_abs_k_ul: float = 0.5
    dysplasia: DysplasiaFeatures = field(default_factory=DysplasiaFeatures)
    iron_store_grade: Optional[IronStoreGrade] = IronStoreGrade.GRADE_3
    cytogenetics_or_mutations: List[str] = field(default_factory=list)
    flow_cytometry_markers: Dict[str, str] = field(default_factory=dict)
    clinical_history: str = ""
    cytopenia_documented: Optional[bool] = None
    persistent_pb_monocytosis_documented: Optional[bool] = None
    plasma_cell_clonality_documented: Optional[bool] = None
    myeloma_defining_event_documented: Optional[bool] = None
    aplastic_anemia_pb_criteria_documented: Optional[bool] = None


@dataclass
class BoneMarrowReport:
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
    who5_interpretation: str = ""
    icc2022_interpretation: str = ""
    interpretation_limitations: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data["cellularity"]["status"] = self.cellularity.status.value
        data["dysplasia_degree"] = self.dysplasia_degree.value
        if self.me_ratio == float("inf"):
            data["me_ratio"] = None
        return data

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, allow_nan=False)


class BoneMarrowDifferentialAnalyzer:
    """Numerical analyzer that returns pattern-level flags rather than autonomous diagnoses."""

    AML_GENETIC_ALIASES = {
        "PML::RARA": ("pml::rara", "pml-rara", "t(15;17)"),
        "RUNX1::RUNX1T1": ("runx1::runx1t1", "runx1-runx1t1", "t(8;21)"),
        "CBFB::MYH11": ("cbfb::myh11", "cbfb-myh11", "inv(16)", "t(16;16)"),
        "KMT2A rearrangement": ("kmt2a", "mll rearrangement", "t(9;11)"),
        "DEK::NUP214": ("dek::nup214", "dek-nup214", "t(6;9)"),
        "MECOM rearrangement": ("mecom", "inv(3)", "t(3;3)"),
        "NUP98 rearrangement": ("nup98",),
        "NPM1 mutation": ("npm1",),
        "CEBPA mutation": ("cebpa",),
        "BCR::ABL1": ("bcr::abl1", "bcr-abl1", "t(9;22)"),
    }

    @staticmethod
    def calculate_me_ratio(granulocytic_count: int, erythroid_count: int) -> float:
        if erythroid_count == 0:
            return float("inf") if granulocytic_count else 0.0
        return round(granulocytic_count / erythroid_count, 2)

    @classmethod
    def _detect_aml_genetics(cls, findings: List[str]) -> List[str]:
        detected = []
        for raw in findings:
            value = raw.strip().lower()
            if any(term in value for term in ("negative", "not detected", "wild type", "wild-type")):
                continue
            for canonical, aliases in cls.AML_GENETIC_ALIASES.items():
                if canonical not in detected and any(alias in value for alias in aliases):
                    detected.append(canonical)
        return detected

    @staticmethod
    def _ipss_r_blast_stratum(blast_pct: float) -> str:
        if blast_pct <= 2:
            return "<= 2% (IPSS-R blast score: 0)"
        if blast_pct < 5:
            return "> 2% to < 5% (IPSS-R blast score: 1)"
        if blast_pct <= 10:
            return "5% to 10% (IPSS-R blast score: 2)"
        return "> 10% (IPSS-R blast score: 3)"

    @classmethod
    def analyze(cls, case: ClinicalCaseInput) -> BoneMarrowReport:
        if not case.case_id.strip():
            raise ValueError("case_id must not be empty")
        if not 0 <= case.peripheral_blood_blast_pct <= 100:
            raise ValueError("peripheral_blood_blast_pct must be between 0 and 100")
        if case.peripheral_blood_monocyte_abs_k_ul < 0:
            raise ValueError("peripheral_blood_monocyte_abs_k_ul must be non-negative")

        count_issues = case.counts.validate()
        dysplasia_issues = case.dysplasia.validate()
        fatal = [
            issue for issue in count_issues + dysplasia_issues
            if issue.startswith(("Count for", "Negative count", "Total nucleated", "erythroid_", "granulocytic_", "megakaryocytic_", "ring_sideroblasts_"))
        ]
        if fatal:
            raise ValueError("; ".join(fatal))
        adequacy_notes = [f"Validation notice: {issue}" for issue in count_issues]

        total = case.counts.total_count()
        blast = case.counts.blast_percentage()
        non_erythroid_blast = case.counts.non_erythroid_blast_percentage()
        me_ratio = case.counts.myeloid_to_erythroid_ratio()
        percentages = case.counts.percentages()
        cellularity = CellularityAssessment(case.patient_age, case.core_cellularity_pct)
        dysplasia_degree = case.dysplasia.get_dysplasia_degree()
        dysplasia_count = case.dysplasia.dysplastic_lineages_count()
        aml_genetics = cls._detect_aml_genetics(case.cytogenetics_or_mutations)

        matched: List[str] = []
        alerts: List[str] = []
        advisories: List[str] = adequacy_notes
        limitations = [
            "This output is an interpretation aid, not a standalone diagnosis or treatment recommendation.",
            "MDS classification requires documented cytopenia plus appropriate exclusions and integrated findings.",
            "WHO-HAEM5 and ICC 2022 differ for several AML/MDS entities and blast thresholds.",
            "The IPSS-R blast stratum is only one component of IPSS-R and applies after MDS is established.",
        ]
        primary = "No major numerical abnormality detected"
        detail = "No disease-level diagnosis assigned from the entered variables"
        who5 = "No WHO-HAEM5 disease entity assigned from the available data."
        icc = "No ICC 2022 disease entity assigned from the available data."

        pb_blast = case.peripheral_blood_blast_pct
        if blast >= 20 or pb_blast >= 20:
            primary = "Acute leukemia blast threshold met"
            detail = "Blasts >=20%; lineage assignment and full acute-leukemia classification are required"
            matched.append(f"Blast threshold >=20% (marrow {blast:.2f}%; peripheral blood {pb_blast:.2f}%).")
            alerts.append("Blasts meet the conventional >=20% acute-leukemia threshold; prompt specialist review is warranted.")
            who5 = ">=20% blasts supports acute-leukemia classification, subject to lineage and entity-specific criteria."
            icc = ">=20% blasts supports AML when myeloid lineage and entity criteria are established."
            if aml_genetics:
                matched.append("AML-associated defining genetic finding(s): " + ", ".join(aml_genetics))
        elif aml_genetics:
            primary = "AML-defining genetic abnormality flag"
            detail = "Apply the exact entity-specific WHO-HAEM5 or ICC blast threshold"
            matched.append("Entered AML-associated genetic finding(s): " + ", ".join(aml_genetics))
            alerts.append("An AML-associated defining genetic abnormality was entered; criteria differ by entity and classification system.")
            who5 = "Several genetically defined AML entities can be classified below the traditional 20% blast threshold; review the exact entity."
            icc = "Most recurrent genetic AML entities use a >=10% blast threshold, with important entity-specific exceptions."
        elif 10 <= blast < 20 or 5 <= pb_blast < 20 or case.dysplasia.auer_rods_present:
            primary = "Increased-blast myeloid neoplasm pattern"
            detail = "WHO-HAEM5 MDS-IB2 range / ICC MDS/AML range only if MDS prerequisites are satisfied"
            matched.append(f"Marrow blasts {blast:.2f}%, PB blasts {pb_blast:.2f}%" + ("; Auer rods entered" if case.dysplasia.auer_rods_present else ""))
            advisories.append("MDS/MDS-AML cannot be assigned from blast count alone; document cytopenia, exclusions and genomic context.")
            who5 = "Values are in the MDS-IB2 range if MDS diagnostic prerequisites and exclusions are satisfied."
            icc = "10-19% blasts are in the MDS/AML range when AML-defining genetics are absent and MDS criteria are otherwise met."
        elif 5 <= blast < 10 or 2 <= pb_blast < 5:
            primary = "Increased-blast marrow pattern"
            detail = "WHO-HAEM5 MDS-IB1 / ICC excess-blast range only if MDS prerequisites are satisfied"
            matched.append(f"Marrow blasts {blast:.2f}%, PB blasts {pb_blast:.2f}%")
            advisories.append("Document cytopenia and exclude mimics before assigning an MDS diagnosis.")
            who5 = "Values are in the MDS-IB1 range if the complete MDS prerequisites are satisfied."
            icc = "Values may fit an MDS excess-blast category if the complete MDS criteria are met."
        elif dysplasia_count:
            primary = "Significant marrow dysplasia pattern"
            detail = f"Dysplasia >=10% in {dysplasia_count} lineage(s); MDS requires additional prerequisites"
            matched.append(
                f"Dysplasia: erythroid {case.dysplasia.erythroid_dysplasia_pct:.1f}%, "
                f"granulocytic {case.dysplasia.granulocytic_dysplasia_pct:.1f}%, "
                f"megakaryocytic {case.dysplasia.megakaryocytic_dysplasia_pct:.1f}%"
            )
            if case.cytopenia_documented is True:
                advisories.append("Cytopenia is documented; correlate with persistence, exclusions, full morphology and genomic findings.")
            else:
                advisories.append("Dysplasia alone is not sufficient for MDS; document cytopenia and exclude secondary causes.")
            who5 = "Morphologic MDS categories require integrated cytopenia, dysplasia, blast assessment and exclusions."
            icc = "MDS classification similarly requires integrated clinical, morphologic and genetic assessment."
        elif percentages.get("plasma_cells", 0) >= 10:
            plasma = percentages["plasma_cells"]
            primary = "Marrow plasmacytosis"
            detail = (
                "Marked plasmacytosis (>=60%); the IMWG biomarker requires clonal plasma cells"
                if plasma >= 60 else
                "Plasma cells >=10%; clonality and myeloma-defining events are required for disease classification"
            )
            matched.append(f"Marrow plasma cells {plasma:.2f}%")
            if case.plasma_cell_clonality_documented and case.myeloma_defining_event_documented:
                primary = "Plasma cell neoplasm criteria documented"
                detail = "Entered context includes clonality and a myeloma-defining event; verify the complete IMWG criteria"
            else:
                advisories.append("Do not diagnose multiple myeloma from the differential percentage alone; establish clonality and CRAB/SLiM criteria.")
            who5 = icc = "Plasma-cell classification requires integrated clonality and clinical criteria beyond the differential count."
        elif cellularity.status == CellularityStatus.SEVERELY_HYPOCELLULAR and blast < 5:
            primary = "Severe marrow hypocellularity"
            if case.aplastic_anemia_pb_criteria_documented is True:
                detail = "Aplastic-anemia-compatible pattern; alternative causes still require exclusion"
                matched.append("Severe hypocellularity plus entered peripheral-blood aplastic-anemia criteria")
            else:
                detail = "Aplastic anemia cannot be assigned without peripheral-blood criteria and exclusion of alternatives"
                matched.append(f"Core cellularity {cellularity.observed_cellularity_pct:.1f}%")
            advisories.append("Correlate with CBC/reticulocytes, exposure history, infection, PNH testing and full marrow review as indicated.")
            who5 = icc = "A marrow-failure pattern is present; aplastic anemia is clinicopathologic and is not established by cellularity alone."

        if case.persistent_pb_monocytosis_documented is True and case.peripheral_blood_monocyte_abs_k_ul >= 0.5:
            advisories.append("Persistent absolute PB monocytosis >=0.5 x10^9/L is entered; CMML assessment also requires >=10% PB monocytes and supporting/clonal and exclusion criteria.")

        if primary == "No major numerical abnormality detected":
            if me_ratio > 4.5:
                primary = "Myeloid-predominant differential pattern"
                detail = f"Elevated M:E ratio ({me_ratio}:1); correlate with morphology and clinical context"
                matched.append(f"M:E ratio {me_ratio}:1")
            elif 0 < me_ratio < 1.2:
                primary = "Erythroid-predominant differential pattern"
                detail = f"Low/inverted M:E ratio ({me_ratio}:1); correlate with erythropoietic response and context"
                matched.append(f"M:E ratio {me_ratio}:1")

        if case.dysplasia.ring_sideroblasts_pct >= 15:
            advisories.append(
                f"Ring sideroblasts {case.dysplasia.ring_sideroblasts_pct:.1f}% entered. "
                "WHO-HAEM5/ICC SF3B1- and ring-sideroblast-associated categories require integrated context."
            )
        elif case.dysplasia.sf3b1_mutation_detected:
            advisories.append("SF3B1 mutation entered; classification requires additional morphologic/genetic prerequisites and exclusions.")
        if case.iron_store_grade in (IronStoreGrade.GRADE_0, IronStoreGrade.GRADE_1):
            advisories.append("Stainable iron is entered as absent/severely decreased; correlate with systemic iron studies.")
        elif case.iron_store_grade in (IronStoreGrade.GRADE_5, IronStoreGrade.GRADE_6):
            advisories.append("Stainable iron is entered as markedly increased; correlate with transfusion history and systemic iron studies.")

        return BoneMarrowReport(
            case.case_id, case.patient_age, total, blast, non_erythroid_blast, me_ratio,
            cellularity, dysplasia_degree, dysplasia_count, primary, detail, matched,
            cls._ipss_r_blast_stratum(blast), alerts, advisories, percentages,
            who5, icc, limitations,
        )


def format_clinical_report(report: BoneMarrowReport) -> str:
    ratio = "inf" if report.me_ratio == float("inf") else f"{report.me_ratio:.2f}"
    lines = [
        "=" * 78,
        f" BONE MARROW DIFFERENTIAL ANALYSIS : {report.case_id}",
        " Interpretation aid — not a standalone diagnosis or treatment recommendation",
        "=" * 78,
        f"Patient age: {report.patient_age} yrs | Total counted: {report.total_cells_counted} cells",
        f"Marrow blasts: {report.marrow_blast_pct:.2f}% | Non-erythroid blasts: {report.non_erythroid_blast_pct:.2f}%",
        f"Myeloid:Erythroid (M:E) ratio: {ratio}:1",
        f"Core cellularity: {report.cellularity.observed_cellularity_pct:.1f}% (approx. age reference {report.cellularity.expected_cellularity_pct:.1f}%) -> {report.cellularity.status.value}",
        f"Dysplasia: {report.dysplasia_degree.value} ({report.dysplastic_lineages_count} lineage(s))",
        f"IPSS-R blast stratum (reference only): {report.ipss_r_blast_score_category}",
        "-" * 78,
        f"INTERPRETIVE CATEGORY: {report.primary_diagnostic_category}",
        f"DETAIL: {report.subclassification}",
        f"WHO-HAEM5: {report.who5_interpretation}",
        f"ICC 2022: {report.icc2022_interpretation}",
        "-" * 78,
    ]
    if report.who_2022_criteria_matched:
        lines += ["Entered findings / thresholds:"] + [f"  - {x}" for x in report.who_2022_criteria_matched]
    if report.critical_alerts:
        lines += ["\nReview flags:"] + [f"  * {x}" for x in report.critical_alerts]
    if report.advisory_recommendations:
        lines += ["\nCorrelation notes:"] + [f"  * {x}" for x in report.advisory_recommendations]
    lines.append("\nDifferential (%):")
    for cell_type, pct in report.differential_percentages.items():
        if pct > 0:
            lines.append(f"  - {cell_type.replace('_', ' ').title():<32}: {pct:>6.2f}%")
    lines += ["\nLimitations:"] + [f"  * {x}" for x in report.interpretation_limitations] + ["=" * 78]
    return "\n".join(lines)
