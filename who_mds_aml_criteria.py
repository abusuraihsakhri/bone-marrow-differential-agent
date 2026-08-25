#!/usr/bin/env python3
"""
WHO 2022 MDS/AML Criteria Validator for Bone Marrow Differential Agent.
Validates morphological findings against WHO 2022 classification criteria
for Myelodysplastic Syndromes and Acute Myeloid Leukemia.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass


@dataclass
class MorphologyFinding:
    """Single morphological finding from bone marrow examination."""
    cell_line: str
    finding: str
    percentage: float = 0.0
    severity: str = "NORMAL"


WHO_MDS_CRITERIA = {
    "MDS_SLIC": {
        "name": "MDS with Single Lineage Dysplasia",
        "criteria": ["dysplasia_in_one_lineage", "blast_pct_lt_5", "ring_sideroblasts_lt_15"],
        "minimum_dysplasia_pct": 10,
    },
    "MDS_ML": {
        "name": "MDS with Multi-Lineage Dysplasia",
        "criteria": ["dysplasia_in_two_plus_lineages", "blast_pct_lt_5", "ring_sideroblasts_lt_15"],
        "minimum_dysplasia_pct": 10,
    },
    "MDS_RS_SL": {
        "name": "MDS with Ring Sideroblasts, Single Lineage",
        "criteria": ["ring_sideroblasts_gte_15", "dysplasia_in_one_lineage", "blast_pct_lt_5"],
    },
    "MDS_RS_ML": {
        "name": "MDS with Ring Sideroblasts, Multi-Lineage",
        "criteria": ["ring_sideroblasts_gte_15", "dysplasia_in_two_plus_lineages", "blast_pct_lt_5"],
    },
    "MDS_EXCESS_BLASTS_1": {
        "name": "MDS with Excess Blasts-1",
        "criteria": ["blast_pct_5_to_9"],
    },
    "MDS_EXCESS_BLASTS_2": {
        "name": "MDS with Excess Blasts-2",
        "criteria": ["blast_pct_10_to_19"],
    },
    "MDS_5Q": {
        "name": "MDS with del(5q)",
        "criteria": ["del_5q_cytogenetics", "blast_pct_lt_5", "normal_aeql_count_or_mild_anemia"],
    },
    "AML_MRC": {
        "name": "AML with Myelodysplasia-Related Changes",
        "criteria": ["dysplasia_in_two_plus_lineages", "blast_pct_gte_20", "no_aml_defining_translocation"],
    },
}


def validate_mds_criteria(findings: List[MorphologyFinding], blast_pct: float,
                          cytogenetics: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Validate bone marrow findings against WHO 2022 MDS criteria."""
    cytogenetics = cytogenetics or {}
    matched_subtypes = []

    dysplastic_lineages = set()
    for finding in findings:
        if finding.percentage >= 10 and finding.finding in ("dysplasia", "dysplastic"):
            dysplastic_lineages.add(finding.cell_line)

    ring_sideroblasts_pct = 0.0
    for finding in findings:
        if "ring_sideroblast" in finding.finding.lower():
            ring_sideroblasts_pct = finding.percentage

    for subtype_id, criteria in WHO_MDS_CRITERIA.items():
        match_count = 0
        total_criteria = len(criteria["criteria"])

        for criterion in criteria["criteria"]:
            if criterion == "dysplasia_in_one_lineage" and len(dysplastic_lineages) >= 1:
                match_count += 1
            elif criterion == "dysplasia_in_two_plus_lineages" and len(dysplastic_lineages) >= 2:
                match_count += 1
            elif criterion == "blast_pct_lt_5" and blast_pct < 5.0:
                match_count += 1
            elif criterion == "blast_pct_5_to_9" and 5.0 <= blast_pct <= 9.0:
                match_count += 1
            elif criterion == "blast_pct_10_to_19" and 10.0 <= blast_pct <= 19.0:
                match_count += 1
            elif criterion == "blast_pct_gte_20" and blast_pct >= 20.0:
                match_count += 1
            elif criterion == "ring_sideroblasts_gte_15" and ring_sideroblasts_pct >= 15.0:
                match_count += 1
            elif criterion == "ring_sideroblasts_lt_15" and ring_sideroblasts_pct < 15.0:
                match_count += 1
            elif criterion == "del_5q_cytogenetics" and cytogenetics.get("del_5q", False):
                match_count += 1
            elif criterion == "no_aml_defining_translocation" and not cytogenetics.get("aml_translocation", False):
                match_count += 1

        if match_count >= total_criteria * 0.8:
            matched_subtypes.append({
                "subtype_id": subtype_id,
                "name": criteria["name"],
                "match_score": match_count / total_criteria,
                "matched_criteria": match_count,
                "total_criteria": total_criteria,
            })

    matched_subtypes.sort(key=lambda x: x["match_score"], reverse=True)

    return {
        "dysplastic_lineages": list(dysplastic_lineages),
        "ring_sideroblasts_pct": ring_sideroblasts_pct,
        "blast_pct": blast_pct,
        "matched_subtypes": matched_subtypes,
        "primary_classification": matched_subtypes[0]["name"] if matched_subtypes else "INDETERMINATE",
        "confidence": matched_subtypes[0]["match_score"] if matched_subtypes else 0.0,
    }


class MdsAmlCriteriaAgent:
    """Sub-agent for WHO 2022 MDS/AML criteria validation."""

    def __init__(self):
        self.agent_name = "MdsAmlCriteriaAgent"

    def evaluate(self, findings: List[MorphologyFinding], blast_pct: float,
                 cytogenetics: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Validate findings against WHO criteria."""
        result = validate_mds_criteria(findings, blast_pct, cytogenetics)
        alerts = []

        if result["primary_classification"] != "INDETERMINATE":
            alerts.append({
                "type": "WHO_CLASSIFICATION_MATCH",
                "severity": "WARNING",
                "message": f"Findings consistent with {result['primary_classification']} "
                           f"(confidence: {result['confidence']:.0%}).",
                "recommendation": "Confirm with cytogenetics and molecular studies."
            })

        if blast_pct >= 20.0:
            alerts.append({
                "type": "AML_BLAST_THRESHOLD",
                "severity": "CRITICAL",
                "message": f"Blast count {blast_pct:.1f}% meets AML diagnostic threshold (>=20%).",
                "recommendation": "Initiate AML workup. Consider flow cytometry and molecular profiling."
            })
        elif blast_pct >= 10.0:
            alerts.append({
                "type": "MDS_EXCESS_BLASTS",
                "severity": "WARNING",
                "message": f"Blast count {blast_pct:.1f}% indicates excess blasts.",
                "recommendation": "Close monitoring. Consider transformation risk assessment."
            })

        return {
            "classification": result,
            "alerts": alerts,
            "dysplastic_lineages": result["dysplastic_lineages"],
            "blast_pct": blast_pct,
        }
