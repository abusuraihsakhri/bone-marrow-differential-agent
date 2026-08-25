#!/usr/bin/env python3
"""
ELN 2022 Risk Stratification for Bone Marrow Differential Agent.
Stratifies AML patients into ELN 2022 risk categories based on
cytogenetics, molecular markers, and blast count.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass


ELN_2022_RISK_CATEGORIES = {
    "favorable": {
        "cytogenetics": [
            "t(8;21)(RUNX1-RUNX1T1)", "inv(16)(p13.1q22)/t(16;16)(CBFB-MYH11)",
            "t(9;11)(KMT2A-MLLT3)", "NPM1_mutation_no_FTLD3",
        ],
        "molecular": ["NPM1_mutation_only", "CBFB-MYH11", "RUNX1-RUNX1T1", "CEBPA_bZIP_in_frame"],
        "description": "Favorable risk with expected >60% 5-year OS",
    },
    "intermediate": {
        "cytogenetics": ["normal_karyotype", "other_abnormalities"],
        "molecular": ["NPM1_wildtype", "FLT3_ITD_low", "FLT3_ITD_high_no_NPM1",
                       "IDH1_R132", "IDH2_R140", "IDH2_R172", "DNMT3A_R882"],
        "description": "Intermediate risk with 40-60% 5-year OS",
    },
    "adverse": {
        "cytogenetics": [
            "complex_karyotype_3_plus", "monosomy_5", "monosomy_7",
            "del_5q", "del_7q", "t(v;11q23)_other", "inv_3q",
            "t(6;9)(DEK-NUP214)", "t(3;3)(GATA2-MECOM)",
        ],
        "molecular": ["TP53_mutation", "RUNX1_mutation", "ASXL1_mutation",
                       "EZH2_mutation", "BCOR_mutation", "SF3B1_mutation",
                       "U2AF1_Q157", "SRSF2_mutation", "STAG2_mutation"],
        "description": "Adverse risk with expected <40% 5-year OS",
    },
}


def stratify_eln2022(cytogenetics: Dict[str, Any], molecular: Optional[Dict[str, Any]] = None,
                     blast_pct: float = 0.0) -> Dict[str, Any]:
    """Stratify AML patient into ELN 2022 risk category."""
    molecular = molecular or {}

    matched = {"favorable": [], "intermediate": [], "adverse": []}

    for category, info in ELN_2022_RISK_CATEGORIES.items():
        for cyto in info["cytogenetics"]:
            if cyto in cytogenetics and cytogenetics[cyto]:
                matched[category].append(f"cyto:{cyto}")

        for mol in info["molecular"]:
            if mol in molecular and molecular[mol]:
                matched[category].append(f"mol:{mol}")

    if matched["adverse"]:
        risk_category = "ADVERSE"
        matched_items = matched["adverse"]
    elif matched["favorable"]:
        risk_category = "FAVORABLE"
        matched_items = matched["favorable"]
    elif matched["intermediate"]:
        risk_category = "INTERMEDIATE"
        matched_items = matched["intermediate"]
    else:
        risk_category = "UNCLASSIFIED"
        matched_items = []

    risk_info = ELN_2022_RISK_CATEGORIES.get(risk_category.lower(), {})

    return {
        "risk_category": risk_category,
        "matched_factors": matched_items,
        "matched_count": len(matched_items),
        "description": risk_info.get("description", "Unable to classify"),
        "blast_pct": blast_pct,
        "favorable_matches": matched["favorable"],
        "intermediate_matches": matched["intermediate"],
        "adverse_matches": matched["adverse"],
    }


class ELNRiskAgent:
    """Sub-agent for ELN 2022 risk stratification."""

    def __init__(self):
        self.agent_name = "ELNRiskAgent"

    def evaluate(self, cytogenetics: Dict[str, Any], molecular: Optional[Dict[str, Any]] = None,
                 blast_pct: float = 0.0) -> Dict[str, Any]:
        """Evaluate ELN 2022 risk stratification."""
        result = stratify_eln2022(cytogenetics, molecular, blast_pct)
        alerts = []

        if result["risk_category"] == "ADVERSE":
            alerts.append({
                "type": "ELN_ADVERSE_RISK",
                "severity": "CRITICAL",
                "message": f"ELN 2022 ADVERSE risk classification with {result['matched_count']} "
                           f"adverse factors identified.",
                "recommendation": "Consider allogeneic stem cell transplant in first remission. "
                                  "Clinical trial enrollment recommended."
            })
        elif result["risk_category"] == "FAVORABLE":
            alerts.append({
                "type": "ELN_FAVORABLE_RISK",
                "severity": "INFO",
                "message": "ELN 2022 FAVORABLE risk classification.",
                "recommendation": "Consolidation chemotherapy without transplant generally appropriate."
            })

        if result["blast_pct"] >= 20.0:
            alerts.append({
                "type": "AML_DIAGNOSTIC_THRESHOLD",
                "severity": "WARNING",
                "message": f"Blast count {result['blast_pct']:.1f}% meets AML diagnostic threshold.",
                "recommendation": "Complete AML workup including flow cytometry."
            })

        return {
            "eln_result": result,
            "alerts": alerts,
            "risk_category": result["risk_category"],
            "matched_factors": result["matched_factors"],
        }
