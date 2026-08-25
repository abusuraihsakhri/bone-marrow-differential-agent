#!/usr/bin/env python3
"""
Dysplasia Scoring System for Bone Marrow Differential Agent.
Quantifies dysplasia severity across cell lineages using morphological criteria
and generates a composite dysplasia score for MDS risk assessment.
"""

from typing import Dict, Any, List
from dataclasses import dataclass


DYSPLASIA_CRITERIA = {
    "erythroid": {
        "features": [
            "nuclear_budding", "nuclear_fragmentation", "multinuclearity",
            "karyorrhexis", "megaloblastoid_change", "vacuolization",
            "internuclear_bridge", "ring_sideroblast",
        ],
        "severity_weights": {
            "absent": 0,
            "mild": 1,
            "moderate": 2,
            "severe": 3,
        },
    },
    "granulocytic": {
        "features": [
            "hypogranular", "pseudo_pelger_huet", "dohle_body_like",
            "hypersegmented_neutrophils", "cytoplasmic_vacuolization",
            "nuclear_hyposegmentation", "abnormal_chromatin_clumping",
        ],
        "severity_weights": {
            "absent": 0,
            "mild": 1,
            "moderate": 2,
            "severe": 3,
        },
    },
    "megakaryocytic": {
        "features": [
            "micromegakaryocyte", "multinucleated_megakaryocyte",
            "hypolobated_megakaryocyte", "separated_nuclear_fragements",
            "abnormal_cloudy_cytoplasm", "large_mononuclear",
        ],
        "severity_weights": {
            "absent": 0,
            "mild": 1,
            "moderate": 2,
            "severe": 3,
        },
    },
}


@dataclass
class DysplasiaAssessment:
    """Dysplasia assessment for a single cell lineage."""
    lineage: str
    features_present: List[str]
    severity: str
    percentage_affected: float
    feature_count: int


def assess_dysplasia(findings: Dict[str, List[str]], percentages: Dict[str, float]) -> Dict[str, Any]:
    """Assess dysplasia severity across all lineages."""
    assessments = []

    for lineage, features in DYSPLASIA_CRITERIA.items():
        present_features = [f for f in findings.get(lineage, []) if f in features["features"]]
        pct = percentages.get(lineage, 0.0)
        feature_count = len(present_features)

        if feature_count == 0:
            severity = "absent"
        elif feature_count <= 2:
            severity = "mild"
        elif feature_count <= 4:
            severity = "moderate"
        else:
            severity = "severe"

        weight = features["severity_weights"][severity]
        weighted_score = weight * (pct / 100.0)

        assessments.append(DysplasiaAssessment(
            lineage=lineage,
            features_present=present_features,
            severity=severity,
            percentage_affected=pct,
            feature_count=feature_count,
        ))

    lineage_scores = {}
    for a in assessments:
        lineage_scores[a.lineage] = {
            "severity": a.severity,
            "feature_count": a.feature_count,
            "percentage_affected": a.percentage_affected,
            "features": a.features_present,
        }

    dysplastic_count = sum(1 for a in assessments if a.severity in ("moderate", "severe"))
    total_features = sum(a.feature_count for a in assessments)

    if dysplastic_count >= 2:
        overall_severity = "multi-lineage"
    elif dysplastic_count == 1:
        overall_severity = "single-lineage"
    else:
        overall_severity = "none"

    composite_score = sum(
        DYSPLASIA_CRITERIA[a.lineage]["severity_weights"][a.severity] * (a.percentage_affected / 100.0)
        for a in assessments
    )

    return {
        "lineage_scores": lineage_scores,
        "dysplastic_lineages": dysplastic_count,
        "overall_severity": overall_severity,
        "composite_dysplasia_score": round(composite_score, 3),
        "total_dysplastic_features": total_features,
    }


class DysplasiaScoringAgent:
    """Sub-agent for dysplasia scoring."""

    def __init__(self):
        self.agent_name = "DysplasiaScoringAgent"

    def evaluate(self, findings: Dict[str, List[str]], percentages: Dict[str, float]) -> Dict[str, Any]:
        """Evaluate dysplasia scoring."""
        result = assess_dysplasia(findings, percentages)
        alerts = []

        if result["overall_severity"] == "multi-lineage":
            alerts.append({
                "type": "MULTI_LINEAGE_DYSPLASIA",
                "severity": "WARNING",
                "message": f"Dysplasia detected in {result['dysplastic_lineages']} lineages "
                           f"(composite score: {result['composite_dysplasia_score']:.3f}).",
                "recommendation": "Consistent with MDS multi-lineage dysplasia. Cytogenetics recommended."
            })
        elif result["overall_severity"] == "single-lineage":
            alerts.append({
                "type": "SINGLE_LINEAGE_DYSPLASIA",
                "severity": "ADVISORY",
                "message": f"Dysplasia in one lineage (score: {result['composite_dysplasia_score']:.3f}).",
                "recommendation": "Monitor for progression. Consider repeat biopsy if clinical suspicion."
            })

        return {
            "dysplasia_result": result,
            "alerts": alerts,
            "severity": result["overall_severity"],
            "composite_score": result["composite_dysplasia_score"],
        }
