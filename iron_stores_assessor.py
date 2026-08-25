#!/usr/bin/env python3
"""
Iron stores assessment with ring sideroblast detection (WHO 2022 / ICCS).

Classifies marrow iron status from Prussian blue grading, ring sideroblast
percentage, SF3B1 mutation status, ferritin and transferrin saturation.

WHO 2022 MDS-RS rule:
    ring sideroblasts >= 15% of erythroid precursors          -> MDS-RS
    ring sideroblasts >= 5% AND SF3B1 mutated                 -> MDS-RS
    (in the presence of another MDS-defining criterion)

Stainable iron interpretation:
    absent            -> iron deficiency (verify with ferritin/TSAT)
    normal/decreased  -> adequate or low stores
    increased         -> iron overload (transfusional siderosis if ferritin high)
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List


@dataclass
class IronAssessment:
    ring_sideroblasts_pct: float = 0.0
    sf3b1_mutated: bool = False
    prussian_blue_grade: str = "normal"   # absent | decreased | normal | increased
    serum_ferritin_ng_ml: float = 100.0
    tsat_pct: float = 25.0
    mds_defining_criteria_present: bool = False

    def __post_init__(self):
        self.prussian_blue_grade = self.prussian_blue_grade.lower()


def assess_iron_stores(a: IronAssessment) -> Dict[str, Any]:
    findings: List[str] = []
    rs = a.ring_sideroblasts_pct

    # --- WHO 2022 ring sideroblast rule -----------------------------------
    if rs >= 15.0 and a.mds_defining_criteria_present:
        classification = "MDS-RS"
        findings.append(f"Ring sideroblasts {rs:.0f}% >= 15% of erythroid precursors")
        molecular = "SF3B1 mutation strongly associated (>80%) with MDS-RS phenotype"
    elif rs >= 5.0 and a.sf3b1_mutated and a.mds_defining_criteria_present:
        classification = "MDS-RS"
        findings.append(f"Ring sideroblasts {rs:.0f}% >= 5% with mutated SF3B1 "
                        "(lowered threshold per WHO 2022)")
        molecular = "SF3B1 mutation supports MDS-RS at reduced RS threshold"
    elif rs >= 15.0:
        classification = "Ring sideroblasts present - MDS-RS pending other criteria"
        findings.append(f"RS {rs:.0f}% meets morphologic threshold; confirm an MDS criterion")
        molecular = "Test SF3B1; RS >= 15% qualifies once any MDS criterion is met"
    elif rs > 0:
        classification = "Normal / non-diagnostic ring sideroblast count"
        findings.append(f"RS {rs:.0f}% below both WHO thresholds (<5% with SF3B1, <15% without)")
        molecular = "N/A"
    else:
        classification = "No ring sideroblasts detected"
        molecular = "N/A"

    # --- Stainable iron ----------------------------------------------------
    if a.prussian_blue_grade == "absent":
        iron_status = "Absent stainable iron - iron deficiency"
        findings.append("Prussian blue negative: no storage iron in macrophages")
        if a.serum_ferritin_ng_ml < 30 or a.tsat_pct < 20:
            findings.append("Ferritin <30 ng/mL or TSAT <20% confirms iron deficiency")
    elif a.prussian_blue_grade == "decreased":
        iron_status = "Decreased iron stores"
        findings.append("Reduced but present storage iron")
    elif a.prussian_blue_grade == "increased":
        if a.serum_ferritin_ng_ml >= 1000:
            iron_status = "Increased iron - transfusional siderosis likely"
            findings.append(f"Ferritin {a.serum_ferritin_ng_ml:.0f} ng/mL with increased "
                            "stainable iron: consider iron chelation review")
        else:
            iron_status = "Increased iron stores"
            findings.append("Elevated macrophage iron without overt overload level")
    else:
        iron_status = "Normal iron stores"
        findings.append("Adequate Prussian blue staining")

    next_steps: List[str] = []
    if "MDS-RS" in classification:
        next_steps += ["Send SF3B1 NGS panel if not done",
                       "Correlate with dysplasia scoring across lineages"]
    if a.prussian_blue_grade == "absent":
        next_steps += ["Iron replacement before ESA therapy",
                       "Evaluate for GI blood loss"]
    if a.serum_ferritin_ng_ml >= 1000:
        next_steps.append("Quantify hepatic/cardiac iron if clinically indicated")

    return {
        "ring_sideroblasts_pct": round(rs, 1),
        "sf3b1_mutated": a.sf3b1_mutated,
        "prussian_blue_grade": a.prussian_blue_grade,
        "iron_status": iron_status,
        "classification": classification,
        "molecular_correlate": molecular,
        "findings": findings,
        "next_steps": next_steps,
    }


if __name__ == "__main__":
    cases = [
        ("MDS-RS, SF3B1+ low RS", IronAssessment(7.0, True, "normal", 180, 28, True)),
        ("Classic MDS-RS", IronAssessment(22.0, False, "normal", 420, 35, True)),
        ("Iron deficiency", IronAssessment(0.0, False, "absent", 8, 9, False)),
        ("Transfusion overload", IronAssessment(2.0, False, "increased", 1650, 55, False)),
        ("Normal", IronAssessment(1.0, False, "normal", 120, 26, False)),
    ]
    for name, case in cases:
        print(f"\n=== {name} ===")
        res = assess_iron_stores(case)
        print(f"classification : {res['classification']}")
        print(f"iron status    : {res['iron_status']}")
        print(f"molecular      : {res['molecular_correlate']}")
