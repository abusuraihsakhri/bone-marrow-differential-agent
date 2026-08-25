#!/usr/bin/env python3
"""
Flow Cytometry Analysis for Bone Marrow Differential Agent.
Analyzes flow cytometry data to identify blast populations,
aberrant immunophenotypes, and minimal residual disease markers.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass


BLAST_MARKERS = ["CD34", "CD117", "CD13", "CD33", "CD19", "CD56", "CD7"]
ABERRANT_PATTERNS = {
    "aberrant_cd7": {"marker": "CD7", "expected_lineage": "T-cell",
                      "context": "AML", "significance": "Aberrant T-antigen on myeloid blasts"},
    "aberrant_cd19": {"marker": "CD19", "expected_lineage": "B-cell",
                       "context": "AML", "significance": "B-lymphoid antigen on AML blasts"},
    "aberrant_cd56": {"marker": "CD56", "expected_lineage": "NK-cell",
                       "context": "AML", "significance": "NK-antigen on AML blasts"},
    "dim_cd34": {"marker": "CD34", "expression": "dim",
                  "context": "AML", "significance": "Dim CD34 suggests aberrant differentiation"},
}


@dataclass
class FlowMarker:
    """Flow cytometry marker measurement."""
    marker: str
    expression: str  # "positive", "negative", "dim", "bright", "heterogeneous"
    percentage: float = 0.0


def analyze_flow_cytometry(markers: List[FlowMarker], cd34_pct: float = 0.0,
                           blast_gate_pct: float = 0.0) -> Dict[str, Any]:
    """Analyze flow cytometry data for blast characterization."""
    marker_dict = {m.marker: m for m in markers}
    blast_markers_present = [m for m in markers if m.marker in BLAST_MARKERS]
    aberrant = []

    for key, pattern in ABERRANT_PATTERNS.items():
        marker = marker_dict.get(pattern["marker"])
        if marker and marker.expression in ("positive", "bright", "dim"):
            aberrant.append({
                "pattern": key,
                "marker": pattern["marker"],
                "significance": pattern["significance"],
            })

    blast_count = cd34_pct if cd34_pct > 0 else blast_gate_pct
    is_aml = blast_count >= 20.0
    is_mds = 5.0 <= blast_count < 20.0

    lineage = "undetermined"
    if marker_dict.get("CD13") and marker_dict["CD13"].expression in ("positive", "bright"):
        lineage = "myeloid"
    elif marker_dict.get("CD19") and marker_dict["CD19"].expression in ("positive", "bright"):
        lineage = "b-lymphoid"
    elif marker_dict.get("CD3") and marker_dict["CD3"].expression in ("positive", "bright"):
        lineage = "t-lymphoid"

    return {
        "blast_count_pct": blast_count,
        "blast_markers": [{"marker": m.marker, "expression": m.expression} for m in blast_markers_present],
        "aberrant_patterns": aberrant,
        "lineage": lineage,
        "diagnosis_suggestion": "AML" if is_aml else "MDS" if is_mds else "Normal/Reactive",
        "immature_population_present": blast_count > 5.0,
    }


def detect_mrd_markers(markers: List[FlowMarker], previous_aberrant: Optional[List[str]] = None) -> Dict[str, Any]:
    """Detect minimal residual disease markers."""
    previous_aberrant = previous_aberrant or []
    mrd_markers = []

    for marker in markers:
        if marker.marker in previous_aberrant and marker.percentage > 0.01:
            mrd_markers.append({
                "marker": marker.marker,
                "percentage": marker.percentage,
                "above_threshold": marker.percentage >= 0.1,
            })

    mrd_positive = len(mrd_markers) > 0

    return {
        "mrd_markers": mrd_markers,
        "mrd_positive": mrd_positive,
        "mrd_percentage": sum(m["percentage"] for m in mrd_markers),
        "sensitivity": "multi-parameter flow cytometry",
    }


class FlowCytometryAgent:
    """Sub-agent for flow cytometry analysis."""

    def __init__(self):
        self.agent_name = "FlowCytometryAgent"

    def evaluate(self, markers: List[FlowMarker], cd34_pct: float = 0.0,
                 blast_gate_pct: float = 0.0) -> Dict[str, Any]:
        """Evaluate flow cytometry data."""
        result = analyze_flow_cytometry(markers, cd34_pct, blast_gate_pct)
        alerts = []

        if result["aberrant_patterns"]:
            alerts.append({
                "type": "ABERRANT_IMMUNOPHENOTYPE",
                "severity": "WARNING",
                "message": f"{len(result['aberrant_patterns'])} aberrant marker pattern(s) detected.",
                "recommendation": "Review for possible lineage infidelity. Correlate with morphology."
            })

        if result["blast_count_pct"] >= 20:
            alerts.append({
                "type": "AML_BLAST_THRESHOLD",
                "severity": "CRITICAL",
                "message": f"Blast population {result['blast_count_pct']:.1f}% by flow cytometry.",
                "recommendation": "Consistent with AML. Proceed with classification workup."
            })

        return {"flow_result": result, "alerts": alerts}
