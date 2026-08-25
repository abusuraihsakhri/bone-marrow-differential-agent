#!/usr/bin/env python3
"""
Pharmacokinetic Dose Optimization for Bone Marrow Differential Agent.
Models drug clearance from blast morphology patterns and recommends
dose adjustments for chemotherapy agents based on pharmacokinetic parameters.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass


@dataclass
class ChemotherapyAgent:
    """Pharmacokinetic parameters for a chemotherapy agent."""
    name: str
    clearance_formula: str  # e.g., "crcl_based", "bmi_based", "renal_adjusted"
    hepatic_metabolism: bool
    renal_elimination: bool
    half_life_hours: float
    therapeutic_index: str  # "narrow", "moderate", "wide"
    organ_impact: str  # "marrow_suppressive", "hepatotoxic", "nephrotoxic"


CHEMO_AGENTS = {
    "cytarabine": ChemotherapyAgent(
        name="Cytarabine", clearance_formula="crcl_based", hepatic_metabolism=False,
        renal_elimination=True, half_life_hours=2.0, therapeutic_index="narrow",
        organ_impact="marrow_suppressive",
    ),
    "daunorubicin": ChemotherapyAgent(
        name="Daunorubicin", clearance_formula="hepatic_adjusted", hepatic_metabolism=True,
        renal_elimination=False, half_life_hours=40.0, therapeutic_index="narrow",
        organ_impact="cardiotoxic",
    ),
    "idarubicin": ChemotherapyAgent(
        name="Idarubicin", clearance_formula="hepatic_adjusted", hepatic_metabolism=True,
        renal_elimination=False, half_life_hours=24.0, therapeutic_index="narrow",
        organ_impact="cardiotoxic",
    ),
    "decitabine": ChemotherapyAgent(
        name="Decitabine", clearance_formula="crcl_based", hepatic_metabolism=False,
        renal_elimination=True, half_life_hours=0.5, therapeutic_index="narrow",
        organ_impact="marrow_suppressive",
    ),
    "venetoclax": ChemotherapyAgent(
        name="Venetoclax", clearance_formula="hepatic_adjusted", hepatic_metabolism=True,
        renal_elimination=False, half_life_hours=26.0, therapeutic_index="moderate",
        organ_impact="marrow_suppressive",
    ),
}


def calculate_dose_adjustment(agent_name: str, baseline_clearance: float,
                               patient_weight: float, blast_pct: float,
                               organ_function: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """Calculate dose adjustment based on patient parameters and blast burden."""
    agent = CHEMO_AGENTS.get(agent_name)
    if not agent:
        return {"error": f"Unknown agent: {agent_name}"}

    organ_function = organ_function or {}
    adjustment_factor = 1.0
    dose_rationale = []

    blast_burden_factor = 1.0
    if blast_pct >= 50:
        blast_burden_factor = 0.75
        dose_rationale.append("High blast burden (>50%): consider dose reduction")
    elif blast_pct >= 20:
        blast_burden_factor = 0.85
        dose_rationale.append("Moderate blast burden (20-50%): consider dose reduction")

    renal_function = organ_function.get("crcl", 120.0)
    if agent.renal_elimination and renal_function < 60:
        renal_factor = max(0.25, renal_function / 120.0)
        adjustment_factor *= renal_factor
        dose_rationale.append(f"Renal impairment (CrCl {renal_function:.0f}): "
                              f"adjustment factor {renal_factor:.2f}")

    hepatic_function = organ_function.get("alt", 40.0)
    if agent.hepatic_metabolism and hepatic_function > 120:
        hepatic_factor = max(0.5, 1.0 - ((hepatic_function - 120) / 200))
        adjustment_factor *= hepatic_factor
        dose_rationale.append(f"Hepatic dysfunction (ALT {hepatic_function:.0f}): "
                              f"adjustment factor {hepatic_factor:.2f}")

    if agent.therapeutic_index == "narrow":
        adjustment_factor *= blast_burden_factor
    else:
        adjustment_factor *= max(blast_burden_factor, 0.8)

    adjusted_clearance = baseline_clearance * adjustment_factor
    recommended_dose_pct = adjustment_factor * 100

    return {
        "agent": agent.name,
        "therapeutic_index": agent.therapeutic_index,
        "organ_impact": agent.organ_impact,
        "adjustment_factor": round(adjustment_factor, 3),
        "recommended_dose_pct": round(recommended_dose_pct, 1),
        "adjusted_clearance": round(adjusted_clearance, 2),
        "dose_rationale": dose_rationale,
        "monitoring_recommendations": get_monitoring(agent_name, blast_pct),
    }


def get_monitoring(agent_name: str, blast_pct: float) -> List[str]:
    """Get monitoring recommendations based on agent and disease state."""
    recs = []
    recs.append("CBC with differential at each cycle")
    if blast_pct >= 20:
        recs.append("Bone marrow biopsy between cycles 1-2 to assess response")
    agent = CHEMO_AGENTS.get(agent_name)
    if agent:
        if agent.organ_impact == "cardiotoxic":
            recs.append("Echocardiogram at baseline and after 2 cycles")
        if agent.renal_elimination:
            recs.append("BUN/creatinine before each cycle")
        if agent.hepatic_metabolism:
            recs.append("Liver function tests weekly during treatment")
    return recs


class PharmacokineticAgent:
    """Sub-agent for pharmacokinetic dose optimization."""

    def __init__(self):
        self.agent_name = "PharmacokineticAgent"

    def evaluate(self, agent_name: str, patient_weight: float, blast_pct: float,
                 organ_function: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
        """Evaluate PK dose optimization."""
        baseline_clearance = 120.0
        result = calculate_dose_adjustment(agent_name, baseline_clearance, patient_weight,
                                           blast_pct, organ_function)
        alerts = []

        if "error" in result:
            alerts.append({
                "type": "INVALID_AGENT",
                "severity": "ERROR",
                "message": result["error"],
                "recommendation": "Select from available agents: " + ", ".join(CHEMO_AGENTS.keys())
            })
        elif result["adjustment_factor"] < 0.7:
            alerts.append({
                "type": "SIGNIFICANT_DOSE_REDUCTION",
                "severity": "WARNING",
                "message": f"Dose reduced to {result['recommended_dose_pct']:.0f}% of standard "
                           f"(factor: {result['adjustment_factor']:.2f}).",
                "recommendation": "Consider alternative agent or supportive care approach."
            })

        return {"pk_result": result, "alerts": alerts}
