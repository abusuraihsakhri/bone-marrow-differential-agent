#!/usr/bin/env python3
"""
Cytogenetics Pattern Recognition for Bone Marrow Differential Agent.
Identifies recurrent cytogenetic abnormalities and correlates with
WHO 2022 classification and prognosis.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass


CYTOGENETIC_PATTERNS = {
    "t(8;21)(q22;q22.1)": {
        "gene_fusion": "RUNX1-RUNX1T1",
        "who_classification": "AML with t(8;21)(RUNX1-RUNX1T1)",
        "eln_category": "favorable",
        "prognosis": "Good",
        "clinical_significance": "Core-binding factor AML. Excellent prognosis with cytarabine-based consolidation.",
    },
    "inv(16)(p13.1q22)": {
        "gene_fusion": "CBFB-MYH11",
        "who_classification": "AML with inv(16)(p13.1q22)/t(16;16)(CBFB-MYH11)",
        "eln_category": "favorable",
        "prognosis": "Good",
        "clinical_significance": "Core-binding factor AML. Good prognosis. May benefit from high-dose cytarabine.",
    },
    "t(9;11)(p21.3;q23.3)": {
        "gene_fusion": "KMT2A-MLLT3",
        "who_classification": "AML with t(9;11)(KMT2A-MLLT3)",
        "eln_category": "favorable",
        "prognosis": "Intermediate-Favorable",
        "clinical_significance": "KMT2A-rearranged AML. Prognosis depends on co-mutations.",
    },
    "t(6;9)(p22.3;q34.1)": {
        "gene_fusion": "DEK-NUP214",
        "who_classification": "AML with t(6;9)(DEK-NUP214)",
        "eln_category": "adverse",
        "prognosis": "Poor",
        "clinical_significance": "Adverse risk. Consider allogeneic transplant.",
    },
    "complex_karyotype": {
        "gene_fusion": "Multiple",
        "who_classification": "AML with myelodysplasia-related changes",
        "eln_category": "adverse",
        "prognosis": "Poor",
        "clinical_significance": ">=3 unrelated abnormalities. High risk. Transplant in CR1.",
    },
    "del(5q)": {
        "gene_fusion": "APC, NPM1 (minimal region)",
        "who_classification": "MDS with del(5q)",
        "eln_category": "adverse",
        "prognosis": "Intermediate",
        "clinical_significance": "5q- syndrome. Lenalidomide responsive if isolated.",
    },
    "monosomy_5": {
        "gene_fusion": "Multiple lost",
        "who_classification": "AML/MDS with monosomy 5",
        "eln_category": "adverse",
        "prognosis": "Poor",
        "clinical_significance": "Adverse risk. Often associated with therapy-related AML.",
    },
    "monosomy_7": {
        "gene_fusion": "Multiple lost",
        "who_classification": "AML/MDS with monosomy 7",
        "eln_category": "adverse",
        "prognosis": "Poor",
        "clinical_significance": "Adverse risk. Consider transplant.",
    },
    "del(7q)": {
        "gene_fusion": "Multiple",
        "who_classification": "AML/MDS with del(7q)",
        "eln_category": "adverse",
        "prognosis": "Intermediate-Poor",
        "clinical_significance": "Adverse risk if unbalanced.",
    },
}


@dataclass
class CytogeneticAbnormality:
    """Detected cytogenetic abnormality."""
    karyotype: str
    abnormality_type: str
    significance: str
    frequency_pct: float


def analyze_cytogenetics(karyotype_string: str) -> Dict[str, Any]:
    """Parse and analyze cytogenetic abnormalities from a karyotype string."""
    detected = []
    karyotype_lower = karyotype_string.lower()

    for pattern, info in CYTOGENETIC_PATTERNS.items():
        pattern_check = pattern.lower().replace("(", "").replace(")", "").replace(";", "")
        karyo_check = karyotype_lower.replace("(", "").replace(")", "").replace(";", "")

        if pattern_check in karyo_check or info["gene_fusion"].lower() in karyotype_lower:
            detected.append(CytogeneticAbnormality(
                karyotype=pattern,
                abnormality_type=info["abnormality_type"],
                significance=info["prognosis"],
                frequency_pct=0.0,
            ))

    has_complex = karyotype_lower.count(">") >= 3 or karyotype_lower.count("add") >= 3
    if has_complex:
        detected.append(CytogeneticAbnormality(
            karyotype="complex_karyotype",
            abnormality_type="structural",
            significance="Poor",
            frequency_pct=0.0,
        ))

    total_abnormalities = len(detected)

    if total_abnormalities == 0:
        eln_category = "unclassified"
        prognosis = "Indeterminate"
    elif any(d.karyotype in CYTOGENETIC_PATTERNS and CYTOGENETIC_PATTERNS[d.karyotype]["eln_category"] == "adverse"
             for d in detected):
        eln_category = "adverse"
        prognosis = "Poor"
    elif any(d.karyotype in CYTOGENETIC_PATTERNS and CYTOGENETIC_PATTERNS[d.karyotype]["eln_category"] == "favorable"
             for d in detected):
        eln_category = "favorable"
        prognosis = "Good"
    else:
        eln_category = "intermediate"
        prognosis = "Intermediate"

    return {
        "karyotype_string": karyotype_string,
        "detected_abnormalities": [
            {"karyotype": d.karyotype, "type": d.abnormality_type, "significance": d.significance}
            for d in detected
        ],
        "total_abnormalities": total_abnormalities,
        "eln_category": eln_category,
        "overall_prognosis": prognosis,
    }


class CytogeneticsAgent:
    """Sub-agent for cytogenetics pattern recognition."""

    def __init__(self):
        self.agent_name = "CytogeneticsAgent"

    def evaluate(self, karyotype_string: str) -> Dict[str, Any]:
        """Evaluate cytogenetics pattern."""
        result = analyze_cytogenetics(karyotype_string)
        alerts = []

        if result["eln_category"] == "adverse":
            alerts.append({
                "type": "ADVERSE_CYTOGENETICS",
                "severity": "CRITICAL",
                "message": f"Adverse cytogenetic risk with {result['total_abnormalities']} "
                           f"abnormalities. Overall prognosis: {result['overall_prognosis']}.",
                "recommendation": "Allogeneic stem cell transplant recommended in first remission."
            })
        elif result["eln_category"] == "favorable":
            alerts.append({
                "type": "FAVORABLE_CYTOGENETICS",
                "severity": "INFO",
                "message": "Favorable cytogenetic profile identified.",
                "recommendation": "Standard consolidation chemotherapy generally appropriate."
            })

        return {"cytogenetics_result": result, "alerts": alerts}
